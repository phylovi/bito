// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"
#include <iostream>
#include <memory>

#include "eigen_sugar.hpp"
#include "numerical_utils.hpp"

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
  ClearTreeCollectionAssociatedState();
  topology_counter_ = tree_collection_.TopologyCounter();
  // Start by adding the rootsplits.
  for (const auto &iter : SBNMaps::RootsplitCounterOf(topology_counter_)) {
    SafeInsert(indexer_, iter.first, index);
    rootsplits_.push_back(iter.first);
    index++;
  }
  // Now add the PCSSs.
  for (const auto &[parent, child_counter] :
       SBNMaps::PCSSCounterOf(topology_counter_)) {
    SafeInsert(parent_to_range_, parent, {index, index + child_counter.size()});
    for (const auto &child_iter : child_counter) {
      const auto &child = child_iter.first;
      SafeInsert(indexer_, parent + child, index);
      SafeInsert(index_to_child_, index, Bitset::ChildSubsplit(parent, child));
      index++;
    }
  }
  sbn_parameters_.resize(index);
  sbn_parameters_.setOnes();
  psp_indexer_ = PSPIndexer(rootsplits_, indexer_);
  taxon_names_ = tree_collection_.TaxonNames();
}

void SBNInstance::CheckSBNMapsAvailable() {
  if (indexer_.empty() || index_to_child_.empty() || parent_to_range_.empty() ||
      rootsplits_.empty() || taxon_names_.empty()) {
    Failwith("Please call ProcessLoadedTrees to prepare your SBN maps.");
  }
}

StringVector SBNInstance::PrettyIndexer() {
  StringVector pretty_representation(indexer_.size());
  for (const auto &[key, idx] : indexer_) {
    if (idx < rootsplits_.size()) {
      pretty_representation[idx] = key.ToString();
    } else {
      pretty_representation[idx] = key.PCSSToString();
    }
  }
  return pretty_representation;
}

void SBNInstance::PrettyPrintIndexer() {
  auto pretty_representation = PrettyIndexer();
  for (size_t i = 0; i < pretty_representation.size(); i++) {
    std::cout << i << "\t" << pretty_representation[i] << std::endl;
  }
}

void SBNInstance::TrainSimpleAverage() {
  auto indexer_representation_counter = SBNMaps::IndexerRepresentationCounterOf(
      indexer_, topology_counter_, sbn_parameters_.size());
  SBNProbability::SimpleAverage(sbn_parameters_, indexer_representation_counter,
                                rootsplits_.size(), parent_to_range_);
}

EigenVectorXd SBNInstance::TrainExpectationMaximization(double alpha, size_t max_iter,
                                                        double score_epsilon) {
  auto indexer_representation_counter = SBNMaps::IndexerRepresentationCounterOf(
      indexer_, topology_counter_, sbn_parameters_.size());
  return SBNProbability::ExpectationMaximization(
      sbn_parameters_, indexer_representation_counter, rootsplits_.size(),
      parent_to_range_, alpha, max_iter, score_epsilon);
}

EigenVectorXd SBNInstance::CalculateSBNProbabilities() {
  EigenVectorXd sbn_parameters_copy = sbn_parameters_;
  SBNProbability::ProbabilityNormalizeParamsInLog(sbn_parameters_copy,
                                                  rootsplits_.size(), parent_to_range_);
  return SBNProbability::ProbabilityOf(sbn_parameters_copy,
                                       MakeIndexerRepresentations());
}

size_t SBNInstance::SampleIndex(std::pair<size_t, size_t> range) const {
  const auto &[start, end] = range;
  Assert(start < end && end <= sbn_parameters_.size(),
         "SampleIndex given an invalid range.");
  // We do not want to overwrite sbn_parameters so we make a copy.
  EigenVectorXd sbn_parameters_subrange = sbn_parameters_.segment(start, end - start);
  NumericalUtils::ProbabilityNormalizeInLog(sbn_parameters_subrange);
  NumericalUtils::Exponentiate(sbn_parameters_subrange);
  std::discrete_distribution<> distribution(sbn_parameters_subrange.begin(),
                                            sbn_parameters_subrange.end());
  // We have to add on range.first because we have taken a slice of the full
  // array, and the sampler treats the beginning of this slice as zero.
  auto result = start + static_cast<size_t>(distribution(random_generator_));
  Assert(result < end, "SampleIndex sampled a value out of range.");
  return result;
}

// This function samples a tree by first sampling the rootsplit, and then
// calling the recursive form of SampleTopology.
Node::NodePtr SBNInstance::SampleTopology(bool rooted) const {
  // Start by sampling a rootsplit.
  size_t rootsplit_index =
      SampleIndex(std::pair<size_t, size_t>(0, rootsplits_.size()));
  const Bitset &rootsplit = rootsplits_.at(rootsplit_index);
  // The addition below turns the rootsplit into a subsplit.
  auto topology = rooted ? SampleTopology(rootsplit + ~rootsplit)
                         : SampleTopology(rootsplit + ~rootsplit)->Deroot();
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

std::vector<IndexerRepresentation> SBNInstance::MakeIndexerRepresentations() const {
  std::vector<IndexerRepresentation> representations;
  representations.reserve(tree_collection_.trees_.size());
  for (const auto &tree : tree_collection_.trees_) {
    representations.push_back(SBNMaps::IndexerRepresentationOf(
        indexer_, tree.Topology(), sbn_parameters_.size()));
  }
  return representations;
}

std::vector<SizeVectorVector> SBNInstance::MakePSPIndexerRepresentations() const {
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

StringSetVector SBNInstance::StringIndexerRepresentationOf(
    IndexerRepresentation indexer_representation) const {
  auto reversed_indexer = StringReversedIndexer();
  StringSetVector string_sets;
  for (const auto &rooted_representation : indexer_representation) {
    StringSet string_set;
    for (const auto index : rooted_representation) {
      SafeInsert(string_set, reversed_indexer[index]);
    }
    string_sets.push_back(std::move(string_set));
  }
  return string_sets;
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

BlockSpecification::ParameterBlockMap SBNInstance::GetPhyloModelParamBlockMap() {
  return GetEngine()->GetPhyloModelBlockSpecification().ParameterBlockMapOf(
      phylo_model_params_);
}

void SBNInstance::MakeEngine(const EngineSpecification &engine_specification,
                             const PhyloModelSpecification &model_specification) {
  CheckSequencesAndTreesLoaded();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
  engine_ =
      std::make_unique<Engine>(engine_specification, model_specification, site_pattern);
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

void SBNInstance::ClearTreeCollectionAssociatedState() {
  sbn_parameters_.resize(0);
  rootsplits_.clear();
  indexer_.clear();
  index_to_child_.clear();
  parent_to_range_.clear();
  topology_counter_.clear();
}

void SBNInstance::PrepareForPhyloLikelihood(
    const PhyloModelSpecification &model_specification, size_t thread_count,
    const std::vector<BeagleFlags> &beagle_flag_vector, const bool use_tip_states,
    std::optional<size_t> tree_count_option) {
  const EngineSpecification engine_specification{thread_count, beagle_flag_vector,
                                                 use_tip_states};
  MakeEngine(engine_specification, model_specification);
  ResizePhyloModelParams(tree_count_option);
}

void SBNInstance::ResizePhyloModelParams(std::optional<size_t> tree_count_option) {
  size_t tree_count =
      tree_count_option ? *tree_count_option : tree_collection_.TreeCount();
  if (tree_count == 0) {
    Failwith(
        "Please add trees to your instance by sampling or loading before "
        "preparing for phylogenetic likelihood calculation.");
  }
  phylo_model_params_.resize(
      tree_count, GetEngine()->GetPhyloModelBlockSpecification().ParameterCount());
}

std::vector<double> SBNInstance::LogLikelihoods() {
  return GetEngine()->LogLikelihoods(tree_collection_, phylo_model_params_, rescaling_);
}

std::vector<std::pair<double, std::vector<double>>> SBNInstance::BranchGradients() {
  return GetEngine()->BranchGradients(tree_collection_, phylo_model_params_,
                                      rescaling_);
}

EigenVectorXd SBNInstance::TopologyGradients(const EigenMatrixXdRef log_f) {

  // Assumption: This function is called from Python side
  // after the trees (both the topology and the branch lengths) are sampled.
  size_t num_trees = tree_collection_.Trees().size();

  EigenVectorXd gradient_vector(sbn_parameters_.size());
  gradient_vector.setZero();

  // Compute log_q_topologies -- remove this once we confirm how this function is going to interact with Python.
//  for (size_t j = 0; j < num_trees; j++) {
//    auto rooted_representation = SBNMaps::RootedIndexerRepresentationOf(indexer_, tree_collection_.GetTree(j).Topology(), OUT_OF_SAMPLE_IDX);
//    log_q_topologies[j] = SBNProbability::SumOf(sbn_parameters_, rooted_representation, 0.);
//  }

  double log_F = NumericalUtils::LogSum(log_f);
  double hat_L = log_F - log(num_trees);
  EigenMatrixXd tilde_w = log_f.array() - log_F;

  // Prepare calculation common to each element of the gradient vector.
  EigenVectorXd temp_gradient_vector(sbn_parameters_.size());
  temp_gradient_vector.setZero();

  double log_norm = NumericalUtils::LogSum(sbn_parameters_.segment(0, rootsplits_.size()));
  temp_gradient_vector.segment(0, rootsplits_.size()) = sbn_parameters_.segment(0, rootsplits_.size()).array() - log_norm;
  //std::cout << temp_gradient_vector << std::endl;
  for (const auto &[_, range] : parent_to_range_) {
    size_t length = (range.second - range.first);
    auto segment = sbn_parameters_.segment(range.first, length);
    log_norm = NumericalUtils::LogSum(segment);
    temp_gradient_vector.segment(range.first, length) = segment.array() - log_norm;
  }
  std::cout << temp_gradient_vector << std::endl;

  // Loop over the topologies to compute the gradient.
  // TODO: Think about ways to make it more efficient.
  for (size_t j = 0; j < num_trees; j++) {
    auto indexer_representation = SBNMaps::IndexerRepresentationOf(indexer_, tree_collection_.GetTree(j).Topology(), OUT_OF_SAMPLE_IDX);
    // TODO: Think about any potential problem that arises from the trees not being independently sampled.
    // Stochastic VI assumes that each of the rooted trees are sampled independently from the variational approximation.
    // Our sampling procedure is 1) sample a rooted tree, 2) unroot it, 3) get unrooted trees by placing a root
    // at every edge. The trees that we obtain as output of step 3) is not independent.
    for (size_t i = 0; i < indexer_representation.size(); i++) {
      double multiplicative_factor = (hat_L - tilde_w(j,i));
      //std::cout << multiplicative_factor << std::endl;
      std::unordered_set<size_t> rooted_indexer_representaiton(indexer_representation.at(i).begin(), indexer_representation.at(i).end());
      std::cout << rooted_indexer_representaiton << std::endl;
      // We need to update all of the SBN parameters for each topology we sampled.
      for (size_t idx = 0; idx < sbn_parameters_.size(); idx++) {
        if (rooted_indexer_representaiton.count(idx) > 0) {
          gradient_vector(idx) += multiplicative_factor * (1-exp(temp_gradient_vector(idx)));
        } else {
          gradient_vector(idx) += multiplicative_factor * -exp(temp_gradient_vector(idx));
        }
        if (idx == 0) {
          std::cout << exp(temp_gradient_vector(idx)) << std::endl;
          std::cout << gradient_vector(0) << std::endl;
        }
      }
    }
  }

  return gradient_vector;
}

// Here we initialize our static random number generator.
std::random_device SBNInstance::random_device_;
std::mt19937 SBNInstance::random_generator_(SBNInstance::random_device_());
