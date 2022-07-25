// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "tp_engine.hpp"
#include "gp_engine.hpp"

TPEngine::TPEngine(GPDAG &dag, size_t site_pattern_count,
                   const std::string &mmap_file_path, bool using_likelihoods,
                   bool using_parsimony)
    : dag_(dag),
      likelihood_pvs_(mmap_file_path, 0, site_pattern_count, 2.0),
      using_likelihoods_(using_likelihoods),
      parsimony_pvs_(mmap_file_path, 0, site_pattern_count, 2.0),
      using_parsimony_(using_parsimony),
      choice_map_(dag) {
  // Initialize node-based data
  GrowNodeData(dag_.NodeCount(), std::nullopt, std::nullopt, true);
  // Initialize edge-based data
  GrowEdgeData(dag_.EdgeCount(), std::nullopt, std::nullopt, true);
  // Initialize scores.
  InitializeChoiceMap();
}

// ** General Scoring

void TPEngine::InitializeChoiceMap() { choice_map_.SelectFirstEdge(); }

Node::NodePtr TPEngine::GetTopTreeTopologyWithEdge(const EdgeId edge_id) const {
  return choice_map_.ExtractTopology(edge_id);
}

double TPEngine::GetScoreOfTopTree(const EdgeId edge_id) const {
  // !ISSUE #440
  Failwith("Currently no implementation.");
  return 0.0;
}

// ** Scoring by Likelihood

void TPEngine::PopulateRootwardPVLikelihoodForNode(const NodeId node_id) {
  // const auto node = dag_.GetDAGNode(node_id);

  // likelihood_pvs_.GetPV(PSVType::PLeft, node_id) = left_likelihood;
  // likelihood_pvs_.GetPV(PSVType::PRight, node_id) = right_likelihood;
}

void TPEngine::PopulateLeafwardPVLikelihoodForNode(const NodeId node_id) {
  // const auto node = dag_.GetDAGNode(node_id);

  // likelihood_pvs_.GetPV(PSVType::PLeft, node_id) = left_likelihood;
  // likelihood_pvs_.GetPV(PSVType::PRight, node_id) = right_likelihood;
}

void TPEngine::InitializeLikelihood() {
  // Rootward Pass (populate P-pvs)
  const auto rootward_node_ids = dag_.RootwardNodeTraversalTrace(true);
  for (const auto node_id : rootward_node_ids) {
    PopulateRootwardPVLikelihoodForNode(node_id);
  }
  // Leafward Pass (populate Q-pvs)
  const auto leafward_node_ids = dag_.LeafwardNodeTraversalTrace(true);
  for (const auto node_id : leafward_node_ids) {
    PopulateLeafwardPVLikelihoodForNode(node_id);
  }

  // !ISSUE #440
  Failwith("Currently no implementation.");
}

void TPEngine::UpdateDAGAfterAddNodePairByLikelihood(const NNIOperation &nni_op) {
  // !ISSUE #440
  Failwith("Currently no implementation.");
}

// ** Scoring by Parsimony

void TPEngine::InitializeParsimony() {
  // !ISSUE #440
  Failwith("Currently no implementation.");
}

void TPEngine::UpdateDAGAfterAddNodePairByParsimony(const NNIOperation &nni_op) {
  // !ISSUE #440
  Failwith("Currently no implementation.");
}

// ** Parameters

void TPEngine::GrowNodeData(const size_t new_node_count,
                            std::optional<const Reindexer> node_reindexer,
                            std::optional<const size_t> explicit_allocation,
                            const bool on_initialization) {
  const size_t old_node_count = GetNodeCount();
  SetNodeCount(new_node_count);
  // Reallocate more space if needed.
  if ((GetPaddedNodeCount() > GetAllocatedNodeCount()) ||
      explicit_allocation.has_value()) {
    SetAllocatedNodeCount(
        size_t(ceil(double(GetPaddedNodeCount()) * resizing_factor_)));
    if (explicit_allocation.has_value()) {
      Assert(explicit_allocation.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      SetAllocatedNodeCount(explicit_allocation.value() + GetSpareNodeCount());
    }
    if (using_likelihoods_) {
      likelihood_pvs_.Resize(new_node_count, GetAllocatedNodeCount());
    }
    if (using_parsimony_) {
      parsimony_pvs_.Resize(new_node_count, GetAllocatedNodeCount());
    }
  }
  // Reindex work space to realign with DAG.
  if (node_reindexer.has_value()) {
    ReindexEdgeData(node_reindexer.value(), old_node_count);
  }
}

void TPEngine::GrowEdgeData(const size_t new_edge_count,
                            std::optional<const Reindexer> edge_reindexer,
                            std::optional<const size_t> explicit_allocation,
                            const bool on_intialization) {
  const size_t old_edge_count = GetEdgeCount();
  SetEdgeCount(new_edge_count);
  // Reallocate more space if needed.
  if ((GetPaddedEdgeCount() > GetAllocatedEdgeCount()) ||
      explicit_allocation.has_value()) {
    SetAllocatedEdgeCount(
        size_t(ceil(double(GetPaddedEdgeCount()) * resizing_factor_)));
    if (explicit_allocation.has_value()) {
      Assert(explicit_allocation.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      SetAllocatedEdgeCount(explicit_allocation.value() + GetSpareEdgeCount());
    }
    branch_lengths_.conservativeResize(GetAllocatedEdgeCount());
  }
  // Resize to fit without deallocating unused memory.
  branch_lengths_.conservativeResize(GetPaddedEdgeCount());
  // Initialize new work space.
  for (size_t i = old_edge_count; i < GetPaddedEdgeCount(); i++) {
    branch_lengths_[i] = default_branch_length_;
  }
  // Reindex work space to realign with DAG.
  if (edge_reindexer.has_value()) {
    ReindexEdgeData(edge_reindexer.value(), old_edge_count);
  }
}

void TPEngine::ReindexNodeData(const Reindexer &node_reindexer,
                               const size_t old_node_count) {
  Assert(node_reindexer.size() == GetNodeCount(),
         "Node Reindexer is the wrong size for GPEngine.");
  Assert(node_reindexer.IsValid(GetNodeCount()), "Node Reindexer is not valid.");

  // Expand node_reindexer into pv_reindexer.
  Reindexer pv_reindexer =
      likelihood_pvs_.BuildPVReindexer(node_reindexer, old_node_count, GetNodeCount());
  // Reindex data vectors
  likelihood_pvs_.Reindex(pv_reindexer);
}

void TPEngine::ReindexEdgeData(const Reindexer &edge_reindexer,
                               const size_t old_edge_count) {
  Assert(edge_reindexer.size() == GetEdgeCount(),
         "Edge Reindexer is the wrong size for GPEngine.");
  Assert(edge_reindexer.IsValid(GetEdgeCount()),
         "Edge Reindexer is not valid for GPEngine size.");
  // Reindex data vectors.
  Reindexer::ReindexInPlace<EigenVectorXd, double>(branch_lengths_, edge_reindexer,
                                                   GetEdgeCount());
}

void TPEngine::GrowSpareNodeData(const size_t new_node_spare_count) {
  if (new_node_spare_count > GetSpareNodeCount()) {
    SetSpareNodeCount(new_node_spare_count);
    GrowNodeData(GetNodeCount());
  }
}

void TPEngine::GrowSpareEdgeData(const size_t new_edge_spare_count) {
  if (new_edge_spare_count > GetSpareEdgeCount()) {
    SetSpareEdgeCount(new_edge_spare_count);
    GrowEdgeData(GetEdgeCount());
  }
}

void TPEngine::InitializeBranchLengthsByTakingFirst(
    const RootedTreeCollection &tree_collection, const BitsetSizeMap &edge_indexer) {
  size_t unique_edge_count = branch_lengths_.size();
  branch_lengths_.setZero();
  EigenVectorXi observed_edge_counts = EigenVectorXi::Zero(unique_edge_count);
  [&observed_edge_counts, this](size_t edge_idx, const RootedTree &tree,
                                const Node *focal_node) {
    branch_lengths_(edge_idx) += tree.BranchLength(focal_node);
    observed_edge_counts(edge_idx)++;
  };
  for (size_t edge_idx = 0; edge_idx < unique_edge_count; edge_idx++) {
    if (observed_edge_counts(edge_idx) == 0) {
      branch_lengths_(edge_idx) = default_branch_length_;
    } else {
      // Normalize the branch length total using the counts to get a mean branch
      // length.
      branch_lengths_(edge_idx) /= static_cast<double>(observed_edge_counts(edge_idx));
    }
  }
}
