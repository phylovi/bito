// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this class is to hold a DAG that we use to build up the operations for
// the generalized pruning operations.

#ifndef SRC_GP_DAG_HPP_
#define SRC_GP_DAG_HPP_

#include "gp_engine.hpp"
#include "quartet_hybrid_request.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "subsplit_dag_node.hpp"
#include "tidy_subsplit_dag.hpp"

class GPDAG : public TidySubsplitDAG {
 public:
  // We store 6 PLVs per subsplit, and index them according to this enum.
  // The notation is as described in the manuscript, but with a slight shift in the
  // position of the tilde. For example P_HAT_TILDE for a subsplit s is
  // \hat{p}(\tilde{s}).
  enum class PLVType { P, P_HAT, P_HAT_TILDE, R_HAT, R, R_TILDE };

  using TidySubsplitDAG::TidySubsplitDAG;

  // Get the index of a PLV of a given type and with a given index.
  static size_t GetPLVIndexStatic(PLVType plv_type, size_t node_count, size_t src_idx);
  size_t GetPLVIndex(PLVType plv_type, size_t src_idx) const;

  // Optimize branch lengths without handling out of date PLVs.
  [[nodiscard]] GPOperationVector ApproximateBranchLengthOptimization(bool use_gradients) const;
  // Schedule branch length, updating PLVs so they are always up to date.
  [[nodiscard]] GPOperationVector BranchLengthOptimization(bool use_gradients);
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

  QuartetHybridRequest QuartetHybridRequestOf(size_t parent_id, bool rotated,
                                              size_t child_id) const;

 private:
  [[nodiscard]] GPOperationVector LeafwardPass(const SizeVector &visit_order) const;
  [[nodiscard]] GPOperationVector RootwardPass(const SizeVector &visit_order) const;

  void AddPhatOperations(const SubsplitDAGNode *node, bool rotated,
                         GPOperationVector &operations) const;
  void AddRhatOperations(const SubsplitDAGNode *node,
                         GPOperationVector &operations) const;
  void OptimizeSBNParametersForASubsplit(const Bitset &subsplit,
                                         GPOperationVector &operations) const;

  GPOperation RUpdateOfRotated(size_t node_id, bool rotated) const;
  // Compute R_HAT(s) = \sum_{t : s < t} P'(s|t) r(t) prior(s|t)
  void UpdateRHat(size_t node_id, GPOperationVector &operations) const;
  void UpdatePHatComputeLikelihood(size_t node_id, size_t child_node_id, bool rotated,
                                   GPOperationVector &operations) const;
  void OptimizeBranchLengthUpdatePHat(size_t node_id, size_t child_node_id,
                                      bool rotated, bool use_gradients,
                                      GPOperationVector &operations) const;
};

#endif  // SRC_GP_DAG_HPP_
