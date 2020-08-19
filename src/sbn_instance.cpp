// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "sbn_instance.hpp"

#include <iostream>
#include <memory>
#include <unordered_set>

#include "eigen_sugar.hpp"
#include "numerical_utils.hpp"

void SBNInstance::PrintStatus() {
  std::cout << "Status for instance '" << name_ << "':\n";
  if (TreeCount()) {
    std::cout << TreeCount() << " unique tree topologies loaded on " << TaxonCount()
              << " leaves.\n";
  } else {
    std::cout << "No trees loaded.\n";
  }
  std::cout << alignment_.Data().size() << " sequences loaded.\n";
}

// ** Building SBN-related items

void SBNInstance::ProcessLoadedTrees() {
  ClearTreeCollectionAssociatedState();
  topology_counter_ = TopologyCounter();
  size_t index;
  std::tie(rootsplits_, indexer_, index_to_child_, parent_to_range_, index) =
      SBNMaps::BuildIndexerBundle(RootsplitCounterOf(topology_counter_),
                                  PCSPCounterOf(topology_counter_));

  sbn_parameters_.resize(index);
  sbn_parameters_.setOnes();
  psp_indexer_ = PSPIndexer(rootsplits_, indexer_);
  taxon_names_ = TaxonNames();
}

void SBNInstance::CheckTopologyCounter() {
  if (TopologyCounter().empty()) {
    Failwith("Please load some trees into your SBN instance.");
  }
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
      pretty_representation[idx] = key.PCSPToString();
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

std::tuple<StringSizeMap, StringSizePairMap> SBNInstance::GetIndexers() const {
  auto str_indexer = StringifyMap(indexer_);
  auto str_parent_to_range = StringifyMap(parent_to_range_);
  std::string rootsplit("rootsplit");
  SafeInsert(str_parent_to_range, rootsplit, {0, rootsplits_.size()});
  return {str_indexer, str_parent_to_range};
}

StringVector SBNInstance::StringReversedIndexer() const {
  std::vector<std::string> reversed_indexer(indexer_.size());
  for (const auto &[key, idx] : indexer_) {
    if (idx < rootsplits_.size()) {
      reversed_indexer[idx] = key.ToString();
    } else {
      reversed_indexer[idx] = key.PCSPToString();
    }
  }
  return reversed_indexer;
}

void SBNInstance::NormalizeSBNParametersInLog(EigenVectorXdRef sbn_parameters) {
  SBNProbability::ProbabilityNormalizeParamsInLog(sbn_parameters, rootsplits_.size(),
                                                  parent_to_range_);
}

// ** Phylogenetic likelihood

Eigen::Ref<EigenMatrixXd> SBNInstance::GetPhyloModelParams() {
  return phylo_model_params_;
}

BlockSpecification::ParameterBlockMap SBNInstance::GetPhyloModelParamBlockMap() {
  return GetEngine()->GetPhyloModelBlockSpecification().ParameterBlockMapOf(
      phylo_model_params_);
}

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
  size_t tree_count = tree_count_option ? *tree_count_option : TreeCount();
  if (tree_count == 0) {
    Failwith(
        "Please add trees to your instance by sampling or loading before "
        "preparing for phylogenetic likelihood calculation.");
  }
  phylo_model_params_.resize(
      tree_count, GetEngine()->GetPhyloModelBlockSpecification().ParameterCount());
}

// ** I/O

void SBNInstance::ReadFastaFile(std::string fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

void SBNInstance::SetAlignment(Alignment &alignment) {
  alignment_ = std::move(alignment);
}

// ** Protected methods

void SBNInstance::MakeEngine(const EngineSpecification &engine_specification,
                             const PhyloModelSpecification &model_specification) {
  CheckSequencesAndTreesLoaded();
  SitePattern site_pattern(alignment_, TagTaxonMap());
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
  auto topology = rooted
                      ? SBNInstance::SampleTopology(rootsplit + ~rootsplit)
                      : SBNInstance::SampleTopology(rootsplit + ~rootsplit)->Deroot();
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

void SBNInstance::ClearTreeCollectionAssociatedState() {
  sbn_parameters_.resize(0);
  rootsplits_.clear();
  indexer_.clear();
  index_to_child_.clear();
  parent_to_range_.clear();
  topology_counter_.clear();
}

void SBNInstance::PushBackRangeForParentIfAvailable(
    const Bitset &parent, SBNInstance::RangeVector &range_vector) {
  if (parent_to_range_.count(parent) > 0) {
    range_vector.push_back(parent_to_range_.at(parent));
  }
}

// Retrieves range of subsplits for each s|t that appears in the tree
// given by rooted_representation.
SBNInstance::RangeVector SBNInstance::GetSubsplitRanges(
    const SizeVector &rooted_representation) {
  RangeVector subsplit_ranges;
  // PROFILE: should we be reserving here?
  subsplit_ranges.emplace_back(0, rootsplits_.size());
  Bitset root = rootsplits_.at(rooted_representation[0]);
  PushBackRangeForParentIfAvailable(root + ~root, subsplit_ranges);
  PushBackRangeForParentIfAvailable(~root + root, subsplit_ranges);
  // Starting at 1 here because we took care of the rootsplit above (the 0th element).
  for (size_t i = 1; i < rooted_representation.size(); i++) {
    Bitset child = index_to_child_.at(rooted_representation[i]);
    PushBackRangeForParentIfAvailable(child, subsplit_ranges);
    PushBackRangeForParentIfAvailable(child.RotateSubsplit(), subsplit_ranges);
  }
  return subsplit_ranges;
}

// This multiplicative factor is the quantity inside the parentheses in eq:nabla in the
// tex.
EigenVectorXd SBNInstance::CalculateMultiplicativeFactors(
    const EigenVectorXdRef log_f) {
  double tree_count = log_f.size();
  double log_F = NumericalUtils::LogSum(log_f);
  double hat_L = log_F - log(tree_count);
  EigenVectorXd tilde_w = log_f.array() - log_F;
  tilde_w = tilde_w.array().exp();
  return hat_L - tilde_w.array();
}

// This multiplicative factor is the quantity inside the parentheses in eq:nablaVIMCO in
// the tex.
EigenVectorXd SBNInstance::CalculateVIMCOMultiplicativeFactors(
    const EigenVectorXdRef log_f) {
  // Use the geometric mean as \hat{f}(\tau^{-j}, \theta^{-j}), in eq:f_hat in
  // the implementation notes.
  size_t tree_count = log_f.size();
  double log_tree_count = log(tree_count);
  double sum_of_log_f = log_f.sum();
  // This has jth entry \hat{f}_{\bm{\phi},{\bm{\psi}}}(\tau^{-j},\bm{\theta}^{-j}),
  // i.e. the log of the geometric mean of each item other than j.
  EigenVectorXd log_geometric_mean = (sum_of_log_f - log_f.array()) / (tree_count - 1);
  EigenVectorXd per_sample_signal(tree_count);
  // This is a vector of entries that when summed become the parenthetical expression
  // in eq:perSampleLearning.
  EigenVectorXd log_f_perturbed = log_f;
  for (size_t j = 0; j < tree_count; j++) {
    log_f_perturbed(j) = log_geometric_mean(j);
    per_sample_signal(j) =
        log_f_perturbed.redux(NumericalUtils::LogAdd) - log_tree_count;
    // Reset the value.
    log_f_perturbed(j) = log_f(j);
  }
  EigenVectorXd multiplicative_factors = CalculateMultiplicativeFactors(log_f);
  multiplicative_factors -= per_sample_signal;
  return multiplicative_factors;
}

// Here we initialize our static random number generator.
std::random_device SBNInstance::random_device_;
std::mt19937 SBNInstance::random_generator_(SBNInstance::random_device_());
