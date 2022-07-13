// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this class is to hold a DAG that we use to build up the operations for
// the generalized pruning operations. Note that rootsplit PCSPs and the DAG root node
// are excluded from operations.
// GPOperationVectors are then consumed and calculated by GPEngine::ProcessOperations().

#pragma once

#include "gp_engine.hpp"
#include "quartet_hybrid_request.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "subsplit_dag_node.hpp"
#include "tidy_subsplit_dag.hpp"
#include "pv_handler.hpp"

using PLVType = PLVHandler::PLVType;

class GPDAG : public TidySubsplitDAG {
 public:
  using TidySubsplitDAG::TidySubsplitDAG;

  // Get the GPEngine index of given PLV type and given node index.
  size_t GetPLVIndex(PLVType plv_type, NodeId node_id) const;

  // ** GPOperations:
  // These methods generate a serial vector of operations, but perform no computation.
  // These operation vectors can then be computed and consumed by the
  // GPEngine::ProcessOperations().

  // Optimize branch lengths without handling out of date PLVs.
  [[nodiscard]] GPOperationVector ApproximateBranchLengthOptimization() const;
  // Schedule branch length, updating PLVs so they are always up to date.
  [[nodiscard]] GPOperationVector BranchLengthOptimization();
  // Compute likelihood values l(s|t) for each child subsplit s by visiting
  // parent subsplit t and generating Likelihood operations for each PCSP s|t.
  // Compute likelihood values l(s) for each rootsplit s by calling
  // MarginalLikelihood().
  [[nodiscard]] GPOperationVector ComputeLikelihoods() const;
  // Do a complete two-pass traversal to correctly populate the PLVs.
  [[nodiscard]] GPOperationVector PopulatePLVs() const;
  // Fill r-PLVs from leaf nodes to the root nodes.
  [[nodiscard]] GPOperationVector LeafwardPass() const;
  // Compute marginal likelihood.
  [[nodiscard]] GPOperationVector MarginalLikelihood() const;
  // Fill p-PLVs from root nodes to the leaf nodes.
  [[nodiscard]] GPOperationVector RootwardPass() const;
  // Optimize SBN parameters.
  [[nodiscard]] GPOperationVector OptimizeSBNParameters() const;
  // Set r-PLVs to zero.
  [[nodiscard]] GPOperationVector SetLeafwardZero() const;
  // Set rhat(s) = stationary for the rootsplits s.
  [[nodiscard]] GPOperationVector SetRhatToStationary() const;
  // Set p-PLVs to zero.
  [[nodiscard]] GPOperationVector SetRootwardZero() const;

  QuartetHybridRequest QuartetHybridRequestOf(NodeId parent_id, bool is_edge_on_left,
                                              NodeId child_id) const;

 private:
  [[nodiscard]] GPOperationVector LeafwardPass(const NodeIdVector &visit_order) const;
  [[nodiscard]] GPOperationVector RootwardPass(const NodeIdVector &visit_order) const;

  void AddPhatOperations(SubsplitDAGNode node, bool is_edge_on_left,
                         GPOperationVector &operations) const;
  void AddRhatOperations(SubsplitDAGNode node, GPOperationVector &operations) const;
  void OptimizeSBNParametersForASubsplit(const Bitset &subsplit,
                                         GPOperationVector &operations) const;

  GPOperation RUpdateOfRotated(NodeId node_id, bool rotated) const;
  // Compute R_HAT(s) = \sum_{t : s < t} P'(s|t) r(t) prior(s|t)
  void UpdateRHat(NodeId node_id, GPOperationVector &operations) const;
  void UpdatePHatComputeLikelihood(NodeId node_id, NodeId child_node_id,
                                   bool is_edge_on_left,
                                   GPOperationVector &operations) const;
  void OptimizeBranchLengthUpdatePHat(NodeId node_id, NodeId child_node_id,
                                      bool is_edge_on_left,
                                      GPOperationVector &operations) const;
};
