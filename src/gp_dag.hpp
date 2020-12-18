// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this class is to hold a DAG that we use to build up the operations for
// the generalized pruning operations.

#ifndef SRC_GP_DAG_HPP_
#define SRC_GP_DAG_HPP_

#include "subsplit_dag_node.hpp"
#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "subsplit_dag.hpp"

class GPDAG : public SubsplitDAG {
 public:

  enum class PLVType { P, P_HAT, P_HAT_TILDE, R_HAT, R, R_TILDE };
  PLVType RPLVType(bool rotated) const {
    return rotated ? PLVType::R_TILDE : PLVType::R;
  }

  using SubsplitDAG::SubsplitDAG;

  // Get the index of a PLV of a given type and with a given index.
  static size_t GetPLVIndexStatic(PLVType plv_type, size_t node_count, size_t src_idx);
  size_t GetPLVIndex(PLVType plv_type, size_t src_idx) const;

  // Schedule branch length optimization.
  [[nodiscard]] GPOperationVector BranchLengthOptimization() const;
  // Compute likelihood values l(s|t) for each child subsplit s by visiting
  // parent subsplit t and generating Likelihood operations for each PCSP s|t.
  // Compute likelihood values l(s) for each rootsplit s by calling
  // MarginalLikelihood().
  [[nodiscard]] GPOperationVector ComputeLikelihoods() const;
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

 private:

  [[nodiscard]] GPOperationVector LeafwardPass(const SizeVector &visit_order) const;
  [[nodiscard]] GPOperationVector RootwardPass(const SizeVector &visit_order) const;

  void AddPhatOperations(const SubsplitDAGNode *node, bool rotated,
                         GPOperationVector &operations) const;
  void AddRhatOperations(const SubsplitDAGNode *node, GPOperationVector &operations) const;
  void OptimizeSBNParametersForASubsplit(const Bitset &subsplit,
                                         GPOperationVector &operations) const;

  void ScheduleBranchLengthOptimization(size_t node_id,
                                        std::unordered_set<size_t> &visited_nodes,
                                        GPOperationVector &operations) const;
  // This function visits and optimizes branches in depth first fashion.
  // It updates p-PLVs and r-PLVs to reflect/propagate the results
  // of branch length optimization from/to other parts of the tree.
  void ScheduleBranchLengthOptimization(bool is_reverse_postorder,
                                        GPOperationVector &operations) const;
  void UpdateRHat(size_t node_id, bool rotated, GPOperationVector &operations) const;
  void UpdatePHatComputeLikelihood(size_t node_id, size_t child_node_id, bool rotated,
                                   GPOperationVector &operations) const;
  void OptimizeBranchLengthUpdatePHat(size_t node_id, size_t child_node_id,
                                      bool rotated,
                                      GPOperationVector &operations) const;
  void UpdateRPLVs(size_t node_id, GPOperationVector &operations) const;
  void OptimizeBranchLengthsUpdatePHatAndPropagateRPLV(
      const SubsplitDAGNode *node, bool rotated, std::unordered_set<size_t> &visited_nodes,
      GPOperationVector &operations) const;
};

#endif  // SRC_GP_DAG_HPP_
