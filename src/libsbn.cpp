// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"

StringPCSSMap StringPCSSMapOf(PCSSDict d) {
  StringPCSSMap d_str;
  for (const auto &iter : d) {
    d_str[iter.first.ToString()] = StringUInt32MapOf(iter.second.Map());
  }
  return d_str;
}

void SBNInstance::FinalizeBeagleInstances() {
  for (const auto &beagle_instance : beagle_instances_) {
    Assert(beagleFinalizeInstance(beagle_instance) == 0,
           "beagleFinalizeInstance gave nonzero return value!");
  }
  beagle_instances_.clear();
  beagle_leaf_count_ = 0;
  beagle_site_count_ = 0;
  }

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
    uint32_t index = 0;
    auto counter = tree_collection_.TopologyCounter();
    // See above for the definitions of these members.
    sbn_parameters_.clear();
    rootsplits_.clear();
    indexer_.clear();
    index_to_child_.clear();
    parent_to_range_.clear();
    // Start by adding the rootsplits.
    for (const auto &iter : RootsplitCounterOf(counter).Map()) {
      SafeInsert(indexer_, iter.first, index);
      rootsplits_.push_back(iter.first);
      index++;
    }
    rootsplit_index_end_ = index;
    // Now add the PCSSs.
    for (const auto &iter : PCSSCounterOf(counter)) {
      const auto &parent = iter.first;
      const auto &child_counter = iter.second;
      SafeInsert(parent_to_range_, parent,
                 {index, index + child_counter.size()});
      for (const auto &child_iter : child_counter.Map()) {
        const auto &child = child_iter.first;
        SafeInsert(indexer_, parent + child, index);
        SafeInsert(index_to_child_, index,
                   Bitset::ChildSubsplit(parent, child));
        index++;
      }
    }
    sbn_parameters_ = std::vector<double>(index, 1.);
  }

  void SBNInstance::CheckSBNMapsAvailable() {
    if (!indexer_.size() || !index_to_child_.size() ||
        !parent_to_range_.size() || !rootsplits_.size()) {
      Failwith("Please call ProcessLoadedTrees to prepare your SBN maps.");
    }
  }

  void SBNInstance::PrintSupports() {
    std::vector<std::string> to_print(indexer_.size());
    for (const auto &iter : indexer_) {
      if (iter.second < rootsplit_index_end_) {
        to_print[iter.second] = iter.first.ToString();
      } else {
        to_print[iter.second] = iter.first.PCSSToString();
      }
    }
    for (size_t i = 0; i < to_print.size(); i++) {
      std::cout << i << "\t" << to_print[i] << std::endl;
    }
  }

  uint32_t SBNInstance::SampleIndex(std::pair<uint32_t, uint32_t> range) const {
    Assert(range.first < range.second && range.second <= sbn_parameters_.size(),
           "SampleIndex given an invalid range.");
    std::discrete_distribution<> distribution(
        sbn_parameters_.begin() + range.first,
        sbn_parameters_.begin() + range.second);
    // We have to add on range.first because we have taken a slice of the full
    // array, and the sampler treats the beginning of this slice as zero.
    auto result =
        range.first + static_cast<uint32_t>(distribution(random_generator_));
    Assert(result < range.second, "SampleIndex sampled a value out of range.");
    return result;
  }

  // This function samples a tree by first sampling the rootsplit, and then
  // calling the recursive form of SampleTopology.
  Node::NodePtr SBNInstance::SampleTopology() const {
    // Start by sampling a rootsplit.
    uint32_t rootsplit_index =
        SampleIndex(std::pair<uint32_t, uint32_t>(0, rootsplit_index_end_));
    const Bitset &rootsplit = rootsplits_.at(rootsplit_index);
    // The addition below turns the rootsplit into a subsplit.
    auto topology = SampleTopology(rootsplit + ~rootsplit)->Deroot();
    topology->Reid();
    return topology;
  }

  // The input to this function is a parent subsplit (of length 2n).
  Node::NodePtr SBNInstance::SampleTopology(
      const Bitset &parent_subsplit) const {
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

  std::vector<IndexerRepresentation> SBNInstance::GetIndexerRepresentations() {
    std::vector<IndexerRepresentation> representations;
    representations.reserve(tree_collection_.trees_.size());
    for (const auto &tree : tree_collection_.trees_) {
      IndexerRepresentationOf(indexer_, tree.Topology());
      representations.push_back(
          IndexerRepresentationOf(indexer_, tree.Topology()));
    }
    return representations;
  }

  // Get the indexer, but reversed and with bitsets appropriately converted to
  // strings.
  StringVector SBNInstance::StringReversedIndexer() {
    std::vector<std::string> reversed_indexer(indexer_.size());
    for (const auto &iter : indexer_) {
      if (iter.second < rootsplit_index_end_) {
        reversed_indexer[iter.second] = iter.first.ToString();
      } else {
        reversed_indexer[iter.second] = iter.first.PCSSToString();
      }
    }
    return reversed_indexer;
  }

  // Turn an IndexerRepresentation into a string representation of the underying
  // bitsets. This is really just so that we can make a test of indexer
  // representations.
  std::pair<StringSet, StringSetVector>
  SBNInstance::StringIndexerRepresentationOf(
      IndexerRepresentation indexer_representation) {
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

  // ** I/O

  std::tuple<StringUInt32Map, StringUInt32PairMap> SBNInstance::GetIndexers() {
    auto indexer_str = StringUInt32MapOf(indexer_);
    StringUInt32PairMap parent_to_range_str;
    for (const auto &iter : parent_to_range_) {
      SafeInsert(parent_to_range_str, iter.first.ToString(), iter.second);
    }
    std::string rootsplit("rootsplit");
    SafeInsert(parent_to_range_str, rootsplit, {0, rootsplit_index_end_});
    return std::tie(indexer_str, parent_to_range_str);
  }

  // This function is really just for testing-- it recomputes from scratch.
  std::pair<StringUInt32Map, StringPCSSMap> SBNInstance::SplitCounters() {
    auto counter = tree_collection_.TopologyCounter();
    return {StringUInt32MapOf(RootsplitCounterOf(counter).Map()),
            StringPCSSMapOf(PCSSCounterOf(counter))};
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
    alignment_.ReadFasta(fname);
  }

  // ** Phylogenetic likelihood

  void SBNInstance::CheckSequencesAndTreesLoaded() {
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

  void SBNInstance::CheckBeagleDimensions() {
    CheckSequencesAndTreesLoaded();
    if (beagle_instances_.size() == 0) {
      Failwith(
          "Call MakeBeagleInstances to make some instances for likelihood "
          "computation.");
    }
    if (alignment_.SequenceCount() != beagle_leaf_count_ ||
        alignment_.Length() != beagle_site_count_) {
      Failwith(
          "Alignment dimensions for current BEAGLE instances do not "
          "match current alignment. Call MakeBeagleInstances again.");
    }
  }

  void SBNInstance::MakeBeagleInstances(int instance_count) {
    // Start by clearing out any existing instances.
    FinalizeBeagleInstances();
    CheckSequencesAndTreesLoaded();
    beagle_leaf_count_ = alignment_.SequenceCount();
    beagle_site_count_ = alignment_.Length();
    SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
    for (auto i = 0; i < instance_count; i++) {
      auto beagle_instance = beagle::CreateInstance(site_pattern);
      beagle::SetJCModel(beagle_instance);
      beagle_instances_.push_back(beagle_instance);
      beagle::PrepareBeagleInstance(beagle_instance, tree_collection_,
                                    site_pattern);
    }
  }

  std::vector<double> SBNInstance::LogLikelihoods() {
    CheckBeagleDimensions();
    return beagle::LogLikelihoods(beagle_instances_, tree_collection_,
                                  rescaling_);
  }

  std::vector<std::pair<double, std::vector<double>>>
  SBNInstance::BranchGradients() {
    return beagle::BranchGradients(beagle_instances_, tree_collection_,
                                   rescaling_);
  }

// Here we initialize our static random number generator.
std::random_device SBNInstance::random_device_;
std::mt19937 SBNInstance::random_generator_(SBNInstance::random_device_());
