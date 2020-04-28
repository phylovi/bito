// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "unrooted_sbn_instance.hpp"
#include <iostream>
#include <memory>
#include <unordered_set>

#include "eigen_sugar.hpp"
#include "numerical_utils.hpp"

// ** Building SBN-related items

void UnrootedSBNInstance::ProcessLoadedTrees() {
  size_t index = 0;
  ClearTreeCollectionAssociatedState();
  topology_counter_ = tree_collection_.TopologyCounter();
  // Start by adding the rootsplits.
  for (const auto &iter : UnrootedSBNMaps::RootsplitCounterOf(topology_counter_)) {
    SafeInsert(indexer_, iter.first, index);
    rootsplits_.push_back(iter.first);
    index++;
  }
  // Now add the PCSSs.
  for (const auto &[parent, child_counter] :
       UnrootedSBNMaps::PCSSCounterOf(topology_counter_)) {
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

void UnrootedSBNInstance::CheckTopologyCounter() {
  if (topology_counter_.empty()) {
    Failwith("Please load some trees into your SBN instance.");
  }
}

void UnrootedSBNInstance::TrainSimpleAverage() {
  CheckTopologyCounter();
  auto indexer_representation_counter = UnrootedSBNMaps::IndexerRepresentationCounterOf(
      indexer_, topology_counter_, sbn_parameters_.size());
  SBNProbability::SimpleAverage(sbn_parameters_, indexer_representation_counter,
                                rootsplits_.size(), parent_to_range_);
}

EigenVectorXd UnrootedSBNInstance::TrainExpectationMaximization(double alpha,
                                                                size_t max_iter,
                                                                double score_epsilon) {
  CheckTopologyCounter();
  auto indexer_representation_counter = UnrootedSBNMaps::IndexerRepresentationCounterOf(
      indexer_, topology_counter_, sbn_parameters_.size());
  return SBNProbability::ExpectationMaximization(
      sbn_parameters_, indexer_representation_counter, rootsplits_.size(),
      parent_to_range_, alpha, max_iter, score_epsilon);
}

EigenVectorXd UnrootedSBNInstance::CalculateSBNProbabilities() {
  EigenVectorXd sbn_parameters_copy = sbn_parameters_;
  SBNProbability::ProbabilityNormalizeParamsInLog(sbn_parameters_copy,
                                                  rootsplits_.size(), parent_to_range_);
  return SBNProbability::ProbabilityOf(sbn_parameters_copy,
                                       MakeIndexerRepresentations());
}

Node::NodePtr UnrootedSBNInstance::SampleTopology() const {
  return SampleTopology(false);
}

void UnrootedSBNInstance::SampleTrees(size_t count) {
  CheckSBNMapsAvailable();
  auto leaf_count = rootsplits_[0].size();
  // 2n-2 because trees are unrooted.
  auto edge_count = 2 * static_cast<int>(leaf_count) - 2;
  tree_collection_.trees_.clear();
  for (size_t i = 0; i < count; i++) {
    std::vector<double> branch_lengths(static_cast<size_t>(edge_count));
    tree_collection_.trees_.emplace_back(
        UnrootedTree(SampleTopology(), std::move(branch_lengths)));
  }
}

std::vector<UnrootedIndexerRepresentation>
UnrootedSBNInstance::MakeIndexerRepresentations() const {
  std::vector<UnrootedIndexerRepresentation> representations;
  representations.reserve(tree_collection_.trees_.size());
  for (const auto &tree : tree_collection_.trees_) {
    representations.push_back(UnrootedSBNMaps::IndexerRepresentationOf(
        indexer_, tree.Topology(), sbn_parameters_.size()));
  }
  return representations;
}

std::vector<SizeVectorVector> UnrootedSBNInstance::MakePSPIndexerRepresentations()
    const {
  std::vector<SizeVectorVector> representations;
  representations.reserve(tree_collection_.trees_.size());
  for (const auto &tree : tree_collection_.trees_) {
    representations.push_back(psp_indexer_.RepresentationOf(tree.Topology()));
  }
  return representations;
}

DoubleVectorVector UnrootedSBNInstance::SplitLengths() const {
  return psp_indexer_.SplitLengths(tree_collection_);
}

// This function is really just for testing-- it recomputes from scratch.
std::pair<StringSizeMap, StringPCSSMap> UnrootedSBNInstance::SplitCounters() const {
  auto counter = tree_collection_.TopologyCounter();
  return {StringifyMap(UnrootedSBNMaps::RootsplitCounterOf(counter).Map()),
          SBNMaps::StringPCSSMapOf(UnrootedSBNMaps::PCSSCounterOf(counter))};
}

// ** I/O

void UnrootedSBNInstance::ReadNewickFile(std::string fname) {
  Driver driver;
  tree_collection_ =
      UnrootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
}

void UnrootedSBNInstance::ReadNexusFile(std::string fname) {
  Driver driver;
  tree_collection_ =
      UnrootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
}

// ** Phylogenetic likelihood

std::vector<double> UnrootedSBNInstance::LogLikelihoods() {
  return GetEngine()->LogLikelihoods(tree_collection_, phylo_model_params_, rescaling_);
}

std::vector<std::pair<double, std::vector<double>>>
UnrootedSBNInstance::BranchGradients() {
  return GetEngine()->BranchGradients(tree_collection_, phylo_model_params_,
                                      rescaling_);
}

void UnrootedSBNInstance::PushBackRangeForParentIfAvailable(
    const Bitset &parent, UnrootedSBNInstance::RangeVector &range_vector) {
  if (parent_to_range_.count(parent) > 0) {
    range_vector.push_back(parent_to_range_.at(parent));
  }
}

// Retrieves range of subsplits for each s|t that appears in the tree
// given by rooted_representation.
UnrootedSBNInstance::RangeVector UnrootedSBNInstance::GetSubsplitRanges(
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

// This gives the gradient of log q at a specific unrooted topology.
// See eq:gradLogQ in the tex, and TopologyGradients for more information about
// normalized_sbn_parameters_in_log.
EigenVectorXd UnrootedSBNInstance::GradientOfLogQ(
    EigenVectorXdRef normalized_sbn_parameters_in_log,
    const UnrootedIndexerRepresentation &indexer_representation) {
  EigenVectorXd grad_log_q = EigenVectorXd::Zero(sbn_parameters_.size());
  double log_q = DOUBLE_NEG_INF;
  for (const auto &rooted_representation : indexer_representation) {
    if (SBNProbability::IsInSBNSupport(rooted_representation, sbn_parameters_.size())) {
      auto subsplit_ranges = GetSubsplitRanges(rooted_representation);
      // Calculate entries in normalized_sbn_parameters_in_log as needed.
      for (const auto &[begin, end] : subsplit_ranges) {
        if (std::isnan(normalized_sbn_parameters_in_log[begin])) {
          // The entry hasn't been filled yet because it's NaN, so fill it.
          auto sbn_parameters_segment = sbn_parameters_.segment(begin, end - begin);
          double log_sum = sbn_parameters_segment.redux(NumericalUtils::LogAdd);
          // We should be extra careful of NaNs when we are using NaN as a sentinel.
          Assert(std::isfinite(log_sum),
                 "GradientOfLogQ encountered non-finite value during calculation.");
          normalized_sbn_parameters_in_log.segment(begin, end - begin) =
              sbn_parameters_segment.array() - log_sum;
        }
      }
      double log_probability_rooted_tree = SBNProbability::SumOf(
          normalized_sbn_parameters_in_log, rooted_representation, 0.0);
      double probability_rooted_tree = exp(log_probability_rooted_tree);
      // We need to look up the subsplits in the tree.
      // Set representation allows fast lookup.
      std::unordered_set<size_t> rooted_representation_as_set(
          rooted_representation.begin(), rooted_representation.end());
      // Now, we actually perform the eq:gradLogQ calculation.
      for (const auto &[begin, end] : subsplit_ranges) {
        for (size_t pcss_idx = begin; pcss_idx < end; pcss_idx++) {
          double indicator_subsplit_in_rooted_tree =
              static_cast<double>(rooted_representation_as_set.count(pcss_idx) > 0);
          grad_log_q[pcss_idx] += probability_rooted_tree *
                                  (indicator_subsplit_in_rooted_tree -
                                   exp(normalized_sbn_parameters_in_log[pcss_idx]));
        }
      }
      log_q = NumericalUtils::LogAdd(log_q, log_probability_rooted_tree);
    }
  }
  grad_log_q.array() *= exp(-log_q);
  return grad_log_q;
}

EigenVectorXd UnrootedSBNInstance::TopologyGradients(const EigenVectorXdRef log_f,
                                                     bool use_vimco) {
  size_t tree_count = tree_collection_.TreeCount();
  EigenVectorXd gradient_vector = EigenVectorXd::Zero(sbn_parameters_.size());
  EigenVectorXd multiplicative_factors =
      use_vimco ? SBNInstance::CalculateVIMCOMultiplicativeFactors(log_f)
                : SBNInstance::CalculateMultiplicativeFactors(log_f);
  // This variable acts as a cache to store normalized SBN parameters in log.
  // Initialization to DOUBLE_NAN indicates that all entries are empty.
  // It is mutated by GradientOfLogQ.
  EigenVectorXd normalized_sbn_parameters_in_log =
      EigenVectorXd::Constant(sbn_parameters_.size(), DOUBLE_NAN);
  for (size_t i = 0; i < tree_count; i++) {
    const auto indexer_representation = UnrootedSBNMaps::IndexerRepresentationOf(
        indexer_, tree_collection_.GetTree(i).Topology(), sbn_parameters_.size());
    // PROFILE: does it matter that we are allocating another sbn_vector_ sized object?
    EigenVectorXd log_grad_q =
        GradientOfLogQ(normalized_sbn_parameters_in_log, indexer_representation);
    log_grad_q.array() *= multiplicative_factors(i);
    gradient_vector += log_grad_q;
  }
  return gradient_vector;
}
