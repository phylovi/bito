// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"
#include <memory>

void SBNInstance::PrintStatus() {
  std::cout << "Status for instance '" << name_ << "':\n";
  if (tree_collection_.TreeCount()) {
    std::cout << TreeCount() << " unique tree topologies loaded on "
              << tree_collection_.TaxonCount() << " leaves.\n";
  } else {
    std::cout << "No trees loaded.\n";
  }
  std::cout << alignment_.Data().size() << " sequences loaded.\n";
}

// ** Building SBN-related items

void SBNInstance::ProcessLoadedTrees() {
  size_t index = 0;
  auto counter = tree_collection_.TopologyCounter();
  // See above for the definitions of these members.
  sbn_parameters_.clear();
  rootsplits_.clear();
  indexer_.clear();
  index_to_child_.clear();
  parent_to_range_.clear();
  // Start by adding the rootsplits.
  for (const auto &iter : SBNMaps::RootsplitCounterOf(counter)) {
    SafeInsert(indexer_, iter.first, index);
    rootsplits_.push_back(iter.first);
    index++;
  }
  // Now add the PCSSs.
  for (const auto &[parent, child_counter] : SBNMaps::PCSSCounterOf(counter)) {
    SafeInsert(parent_to_range_, parent, {index, index + child_counter.size()});
    for (const auto &child_iter : child_counter) {
      const auto &child = child_iter.first;
      SafeInsert(indexer_, parent + child, index);
      SafeInsert(index_to_child_, index, Bitset::ChildSubsplit(parent, child));
      index++;
    }
  }
  sbn_parameters_ = std::vector<double>(index, 1.);
  psp_indexer_ = PSPIndexer(rootsplits_, indexer_);
}

void SBNInstance::CheckSBNMapsAvailable() {
  if (!indexer_.size() || !index_to_child_.size() || !parent_to_range_.size() ||
      !rootsplits_.size()) {
    Failwith("Please call ProcessLoadedTrees to prepare your SBN maps.");
  }
}

void SBNInstance::PrintSupports() {
  std::vector<std::string> to_print(indexer_.size());
  for (const auto &[key, idx] : indexer_) {
    if (idx < rootsplits_.size()) {
      to_print[idx] = key.ToString();
    } else {
      to_print[idx] = key.PCSSToString();
    }
  }
  for (size_t i = 0; i < to_print.size(); i++) {
    std::cout << i << "\t" << to_print[i] << std::endl;
  }
}

size_t SBNInstance::SampleIndex(std::pair<size_t, size_t> range) const {
  Assert(range.first < range.second && range.second <= sbn_parameters_.size(),
         "SampleIndex given an invalid range.");
  std::discrete_distribution<> distribution(
      // Lordy, these integer types.
      sbn_parameters_.begin() + static_cast<ptrdiff_t>(range.first),
      sbn_parameters_.begin() + static_cast<ptrdiff_t>(range.second));
  // We have to add on range.first because we have taken a slice of the full
  // array, and the sampler treats the beginning of this slice as zero.
  auto result =
      range.first + static_cast<size_t>(distribution(random_generator_));
  Assert(result < range.second, "SampleIndex sampled a value out of range.");
  return result;
}

// This function samples a tree by first sampling the rootsplit, and then
// calling the recursive form of SampleTopology.
Node::NodePtr SBNInstance::SampleTopology() const {
  // Start by sampling a rootsplit.
  size_t rootsplit_index =
      SampleIndex(std::pair<size_t, size_t>(0, rootsplits_.size()));
  const Bitset &rootsplit = rootsplits_.at(rootsplit_index);
  // The addition below turns the rootsplit into a subsplit.
  auto topology = SampleTopology(rootsplit + ~rootsplit)->Deroot();
  topology->Polish();
  return topology;
}

// The input to this function is a parent subsplit (of length 2n).
Node::NodePtr SBNInstance::SampleTopology(const Bitset &parent_subsplit) const {
  auto process_subsplit = [this](const Bitset &parent) {
    auto singleton_option = parent.SplitChunk(1).SingletonOption();
    if (singleton_option) {
      return Node::Leaf(*singleton_option);
    }  // else
    auto child_index = SampleIndex(parent_to_range_.at(parent));
    return SampleTopology(index_to_child_.at(child_index));
  };
  return Node::Join(process_subsplit(parent_subsplit),
                    process_subsplit(parent_subsplit.RotateSubsplit()));
}

void SBNInstance::SampleTrees(size_t count) {
  CheckSBNMapsAvailable();
  auto leaf_count = rootsplits_[0].size();
  // 2n-2 because trees are unrooted.
  auto edge_count = 2 * static_cast<int>(leaf_count) - 2;
  tree_collection_.trees_.clear();
  for (size_t i = 0; i < count; i++) {
    std::vector<double> branch_lengths(static_cast<size_t>(edge_count));
    tree_collection_.trees_.emplace_back(
        Tree(SampleTopology(), std::move(branch_lengths)));
  }
}

std::vector<IndexerRepresentation> SBNInstance::GetIndexerRepresentations()
    const {
  std::vector<IndexerRepresentation> representations;
  representations.reserve(tree_collection_.trees_.size());
  for (const auto &tree : tree_collection_.trees_) {
    representations.push_back(
        SBNMaps::IndexerRepresentationOf(indexer_, tree.Topology()));
  }
  return representations;
}

std::vector<SizeVectorVector> SBNInstance::GetPSPIndexerRepresentations()
    const {
  std::vector<SizeVectorVector> representations;
  representations.reserve(tree_collection_.trees_.size());
  for (const auto &tree : tree_collection_.trees_) {
    representations.push_back(psp_indexer_.RepresentationOf(tree.Topology()));
  }
  return representations;
}

StringVector SBNInstance::StringReversedIndexer() const {
  std::vector<std::string> reversed_indexer(indexer_.size());
  for (const auto &[key, idx] : indexer_) {
    if (idx < rootsplits_.size()) {
      reversed_indexer[idx] = key.ToString();
    } else {
      reversed_indexer[idx] = key.PCSSToString();
    }
  }
  return reversed_indexer;
}

std::pair<StringSet, StringSetVector>
SBNInstance::StringIndexerRepresentationOf(
    IndexerRepresentation indexer_representation) const {
  auto reversed_indexer = StringReversedIndexer();
  auto rootsplit_indices = indexer_representation.first;
  auto pcss_index_vector = indexer_representation.second;
  StringSet rootsplit_string_set;
  for (auto index : rootsplit_indices) {
    SafeInsert(rootsplit_string_set, reversed_indexer[index]);
  }
  StringSetVector pcss_string_sets;
  for (const auto &pcss_indices : pcss_index_vector) {
    StringSet pcss_string_set;
    for (auto index : pcss_indices) {
      SafeInsert(pcss_string_set, reversed_indexer[index]);
    }
    pcss_string_sets.push_back(std::move(pcss_string_set));
  }
  return {rootsplit_string_set, pcss_string_sets};
}

DoubleVectorVector SBNInstance::SplitLengths() const {
  return psp_indexer_.SplitLengths(tree_collection_);
}

// ** I/O

std::tuple<StringSizeMap, StringSizePairMap> SBNInstance::GetIndexers() const {
  auto str_indexer = StringifyMap(indexer_);
  auto str_parent_to_range = StringifyMap(parent_to_range_);
  std::string rootsplit("rootsplit");
  SafeInsert(str_parent_to_range, rootsplit, {0, rootsplits_.size()});
  return {str_indexer, str_parent_to_range};
}

// This function is really just for testing-- it recomputes from scratch.
std::pair<StringSizeMap, StringPCSSMap> SBNInstance::SplitCounters() const {
  auto counter = tree_collection_.TopologyCounter();
  return {StringifyMap(SBNMaps::RootsplitCounterOf(counter).Map()),
          SBNMaps::StringPCSSMapOf(SBNMaps::PCSSCounterOf(counter))};
}

void SBNInstance::ReadNewickFile(std::string fname) {
  Driver driver;
  tree_collection_ = driver.ParseNewickFile(fname);
}

void SBNInstance::ReadNexusFile(std::string fname) {
  Driver driver;
  tree_collection_ = driver.ParseNexusFile(fname);
}

void SBNInstance::ReadFastaFile(std::string fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

// ** Phylogenetic likelihood

void SBNInstance::CheckSequencesAndTreesLoaded() const {
  if (alignment_.SequenceCount() == 0) {
    Failwith(
        "Load an alignment into your SBNInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
  if (TreeCount() == 0) {
    Failwith(
        "Load some trees into your SBNInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
}

Eigen::Ref<EigenMatrixXd> SBNInstance::GetPhyloModelParams() {
  return phylo_model_params_;
}

BlockSpecification::ParameterBlockMap
SBNInstance::GetPhyloModelParamBlockMap() {
  return GetEngine()->GetBlockSpecification().ParameterBlockMapOf(
      phylo_model_params_);
}

void SBNInstance::MakeEngine(PhyloModelSpecification specification,
                             size_t thread_count) {
  CheckSequencesAndTreesLoaded();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
  engine_ = std::make_unique<Engine>(specification, site_pattern, thread_count);
}

Engine *SBNInstance::GetEngine() const {
  if (engine_ != nullptr) {
    return engine_.get();
  }
  // else
  Failwith(
      "Engine not available. Call PrepareForPhyloLikelihood to make an "
      "engine for phylogenetic likelihood computation computation.");
}

void SBNInstance::PrepareForPhyloLikelihood(
    PhyloModelSpecification specification, size_t thread_count,
    std::optional<size_t> tree_count_option) {
  MakeEngine(specification, thread_count);
  ResizePhyloModelParams(tree_count_option);
}

void SBNInstance::ResizePhyloModelParams(
    std::optional<size_t> tree_count_option) {
  size_t tree_count =
      tree_count_option ? *tree_count_option : tree_collection_.TreeCount();
  if (tree_count == 0) {
    Failwith(
        "Please add trees to your instance by sampling or loading before "
        "preparing for phylogenetic likelihood calculation.");
  }
  phylo_model_params_.resize(
      tree_count, GetEngine()->GetBlockSpecification().ParameterCount());
}

std::vector<double> SBNInstance::LogLikelihoods() {
  return GetEngine()->LogLikelihoods(tree_collection_, phylo_model_params_);
}

std::vector<std::pair<double, std::vector<double>>>
SBNInstance::BranchGradients() {
  return GetEngine()->BranchGradients(tree_collection_, phylo_model_params_);
}

// Here we initialize our static random number generator.
std::random_device SBNInstance::random_device_;
std::mt19937 SBNInstance::random_generator_(SBNInstance::random_device_());
