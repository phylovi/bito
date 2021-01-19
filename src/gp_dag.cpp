// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_dag.hpp"
#include "numerical_utils.hpp"

using namespace GPOperations;  // NOLINT

size_t GPDAG::GetPLVIndexStatic(PLVType plv_type, size_t node_count, size_t src_idx) {
  switch (plv_type) {
    case PLVType::P:
      return src_idx;
    case PLVType::P_HAT:
      return node_count + src_idx;
    case PLVType::P_HAT_TILDE:
      return 2 * node_count + src_idx;
    case PLVType::R_HAT:
      return 3 * node_count + src_idx;
    case PLVType::R:
      return 4 * node_count + src_idx;
    case PLVType::R_TILDE:
      return 5 * node_count + src_idx;
    default:
      Failwith("Invalid PLV index requested.");
  }
}

size_t GPDAG::GetPLVIndex(PLVType plv_type, size_t src_idx) const {
  return GetPLVIndexStatic(plv_type, dag_nodes_.size(), src_idx);
}

void GPDAG::UpdateRPLVs(size_t node_id, GPOperationVector &operations) const {
  // Compute R_HAT(s) = \sum_{t : s < t} q(s|t) P(s|t) r(t).
  operations.push_back(Zero{GetPLVIndex(PLVType::R_HAT, node_id)});
  UpdateRHat(node_id, false, operations);
  UpdateRHat(node_id, true, operations);

  // Update r(s) and r_tilde(s)
  operations.push_back(Multiply{GetPLVIndex(PLVType::R, node_id),
                                GetPLVIndex(PLVType::R_HAT, node_id),
                                GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
  operations.push_back(Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                                GetPLVIndex(PLVType::R_HAT, node_id),
                                GetPLVIndex(PLVType::P_HAT, node_id)});
}

GPOperationVector GPDAG::BranchLengthOptimization() const {
  GPOperationVector operations;
  DepthFirstWithAction(SubsplitDAGTraversalAction(
      // BeforeNode
      [this, &operations](size_t node_id) {
        if (!GetDagNode(node_id)->IsRoot()) {
          UpdateRPLVs(node_id, operations);
        }
      },
      // AfterNode
      [this, &operations](size_t node_id) {
        operations.push_back(Multiply{GetPLVIndex(PLVType::P, node_id),
                                      GetPLVIndex(PLVType::P_HAT, node_id),
                                      GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
      },
      // BeforeNodeClade
      [this, &operations](size_t node_id, bool rotated) {
        const PLVType p_hat_plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
        operations.push_back(Zero{GetPLVIndex(p_hat_plv_type, node_id)});
      },
      // AfterNodeClade
      [this, &operations](size_t node_id, bool rotated) {
        const PLVType p_hat_plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
        const PLVType r_plv_type = rotated ? PLVType::R : PLVType::R_TILDE;
        // Depending on rotated, either:
        // r_tilde(t) = r_hat(t) \circ p_hat(t)
        // r(t) = r_hat(t) \circ p_hat_tilde(t)
        operations.push_back(Multiply{GetPLVIndex(r_plv_type, node_id),
                                      GetPLVIndex(PLVType::R_HAT, node_id),
                                      GetPLVIndex(p_hat_plv_type, node_id)});
      },
      // VisitEdge
      [this, &operations](size_t node_id, size_t child_id, bool rotated) {
        OptimizeBranchLengthUpdatePHat(node_id, child_id, rotated, operations);
      }));
  return operations;
}

GPOperationVector GPDAG::OldBranchLengthOptimization() const {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  IterateOverRootsplitIds([this, &operations, &visited_nodes](size_t rootsplit_id) {
    DeprecatedScheduleBranchLengthOptimization(rootsplit_id, visited_nodes, operations);
  });

  return operations;
}

void GPDAG::DeprecatedScheduleBranchLengthOptimization(
    size_t node_id, std::unordered_set<size_t> &visited_nodes,
    GPOperationVector &operations) const {
  visited_nodes.insert(node_id);
  const auto node = GetDagNode(node_id);

  if (node->IsLeaf()) {
    return;
  }

  if (!node->IsRoot()) {
    UpdateRPLVs(node_id, operations);
  }

  DeprecatedOptimizeBranchLengthsUpdatePHatAndPropagateRPLV(node, false, visited_nodes,
                                                            operations);
  DeprecatedOptimizeBranchLengthsUpdatePHatAndPropagateRPLV(node, true, visited_nodes,
                                                            operations);

  // Update p(t).
  operations.push_back(Multiply{GetPLVIndex(PLVType::P, node_id),
                                GetPLVIndex(PLVType::P_HAT, node_id),
                                GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
};

void GPDAG::DeprecatedOptimizeBranchLengthsUpdatePHatAndPropagateRPLV(
    const SubsplitDAGNode *node, bool rotated,
    std::unordered_set<size_t> &visited_nodes, GPOperationVector &operations) const {
  PLVType p_hat_plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
  PLVType r_plv_type = rotated ? PLVType::R : PLVType::R_TILDE;

  size_t node_id = node->Id();

  operations.push_back(Zero{GetPLVIndex(p_hat_plv_type, node_id)});
  IterateOverLeafwardEdges(node, rotated,
                           [this, &operations, &rotated, &visited_nodes,
                            &node_id](const SubsplitDAGNode *child_node) {
                             size_t child_id = child_node->Id();
                             if (visited_nodes.count(child_id) == 0) {
                               DeprecatedScheduleBranchLengthOptimization(
                                   child_id, visited_nodes, operations);
                             }
                             OptimizeBranchLengthUpdatePHat(node_id, child_id, rotated,
                                                            operations);
                           });
  // Depending on rotated, either:
  // r_tilde(t) = r_hat(t) \circ p_hat(t)
  // r(t) = r_hat(t) \circ p_hat_tilde(t)
  operations.push_back(Multiply{GetPLVIndex(r_plv_type, node_id),
                                GetPLVIndex(PLVType::R_HAT, node_id),
                                GetPLVIndex(p_hat_plv_type, node_id)});
}

GPOperationVector GPDAG::ComputeLikelihoods() const {
  GPOperationVector operations;
  IterateOverRealNodes([this, &operations](const SubsplitDAGNode *node) {
    IterateOverLeafwardEdges(
        node, [this, node, &operations](const bool rotated,
                                        const SubsplitDAGNode *child_node) {
          const auto gpcsp_idx =
              GetGPCSPIndex(node->GetBitset(rotated), child_node->GetBitset());
          operations.push_back(Likelihood{gpcsp_idx,
                                          GetPLVIndex(RPLVType(rotated), node->Id()),
                                          GetPLVIndex(PLVType::P, child_node->Id())});
        });
  });

  const auto marginal_likelihood_operations = MarginalLikelihood();
  operations.insert(operations.end(), marginal_likelihood_operations.begin(),
                    marginal_likelihood_operations.end());
  return operations;
}

GPOperationVector GPDAG::LeafwardPass() const {
  return LeafwardPass(LeafwardPassTraversal());
}

GPOperationVector GPDAG::MarginalLikelihood() const {
  GPOperationVector operations;
  for (size_t rootsplit_idx = 0; rootsplit_idx < rootsplits_.size(); rootsplit_idx++) {
    const auto rootsplit = rootsplits_[rootsplit_idx];
    const auto root_subsplit = rootsplit + ~rootsplit;
    size_t root_idx = subsplit_to_id_.at(root_subsplit);
    operations.push_back(GPOperations::IncrementMarginalLikelihood{
        GetPLVIndex(PLVType::R_HAT, root_idx), rootsplit_idx,
        GetPLVIndex(PLVType::P, root_idx)});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass() const {
  return RootwardPass(RootwardPassTraversal());
}

GPOperationVector GPDAG::OptimizeSBNParameters() const {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (size_t &id : LeafwardPassTraversal()) {
    const auto node = GetDagNode(id);
    OptimizeSBNParametersForASubsplit(node->GetBitset(), operations);
    OptimizeSBNParametersForASubsplit(node->GetBitset().RotateSubsplit(), operations);
  }
  operations.push_back(UpdateSBNProbabilities{0, rootsplits_.size()});
  return operations;
}

GPOperationVector GPDAG::SetLeafwardZero() const {
  GPOperationVector operations;
  const auto node_count = dag_nodes_.size();
  for (size_t i = 0; i < node_count; i++) {
    operations.push_back(Zero{GetPLVIndex(PLVType::R_HAT, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::R, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::R_TILDE, i)});
  }
  return operations;
}

GPOperationVector GPDAG::SetRhatToStationary() const {
  GPOperationVector operations;
  IterateOverRootsplitIds([this, &operations](size_t rootsplit_id) {
    size_t root_gpcsp_idx = gpcsp_indexer_.at(GetDagNode(rootsplit_id)->GetBitset());
    operations.push_back(SetToStationaryDistribution{
        GetPLVIndex(PLVType::R_HAT, rootsplit_id), root_gpcsp_idx});
  });
  return operations;
}

GPOperationVector GPDAG::SetRootwardZero() const {
  GPOperationVector operations;
  const auto node_count = dag_nodes_.size();
  for (size_t i = taxon_count_; i < node_count; i++) {
    operations.push_back(Zero{GetPLVIndex(PLVType::P, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT_TILDE, i)});
  }
  return operations;
}

GPOperationVector GPDAG::LeafwardPass(const SizeVector &visit_order) const {
  GPOperationVector operations;
  for (const size_t node_id : visit_order) {
    const auto node = GetDagNode(node_id);

    // Build rhat(s) via rhat(s) += \sum_t q(s|t) P'(s|t) r(t)
    AddRhatOperations(node, operations);
    // Multiply to get r(s) = rhat(s) \circ phat(s_tilde).
    operations.push_back(Multiply{GetPLVIndex(PLVType::R, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
    // Multiply to get r(s_tilde) = rhat(s) \circ phat(s).
    operations.push_back(Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT, node_id)});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass(const SizeVector &visit_order) const {
  GPOperationVector operations;
  for (const size_t node_id : visit_order) {
    const auto node = GetDagNode(node_id);
    if (node->IsLeaf()) {
      continue;
    }
    // Build phat(s).
    AddPhatOperations(node, false, operations);
    // Build phat(s_tilde).
    AddPhatOperations(node, true, operations);
    // Multiply to get p(s) = phat(s) \circ phat(s_tilde).
    operations.push_back(Multiply{node_id, GetPLVIndex(PLVType::P_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
  }
  return operations;
}

GPOperationVector GPDAG::PopulatePLVs() {
  GPOperationVector operations;
  auto append_operations = [&operations](GPOperationVector &&new_operations) {
    std::move(new_operations.begin(), new_operations.end(),
              std::back_inserter(operations));
  };
  append_operations(SetRootwardZero());
  append_operations(SetLeafwardZero());
  append_operations(SetRhatToStationary());
  append_operations(RootwardPass());
  append_operations(LeafwardPass());
  return operations;
}

// Take in some new operations, determine an appropriate PrepForMarginalization for
// them, then append the PrepForMarginalization and the new operations to `operations`
// (in that order).
void AppendOperationsAfterPrepForMarginalization(
    GPOperationVector &operations, const GPOperationVector &new_operations) {
  if (!new_operations.empty()) {
    operations.push_back(PrepForMarginalizationOfOperations(new_operations));
    operations.insert(operations.end(), new_operations.begin(), new_operations.end());
  }
}

void GPDAG::AddPhatOperations(const SubsplitDAGNode *node, bool rotated,
                              GPOperationVector &operations) const {
  PLVType plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
  const auto parent_subsplit = node->GetBitset(rotated);
  const size_t dest_idx = GetPLVIndex(plv_type, node->Id());
  GPOperationVector new_operations;
  for (const size_t &child_idx : node->GetLeafward(rotated)) {
    const auto gpcsp_idx =
        GetGPCSPIndex(parent_subsplit, GetDagNode(child_idx)->GetBitset());
    new_operations.push_back(IncrementWithWeightedEvolvedPLV{
        dest_idx, gpcsp_idx, GetPLVIndex(PLVType::P, child_idx)});
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::AddRhatOperations(const SubsplitDAGNode *node,
                              GPOperationVector &operations) const {
  const auto subsplit = node->GetBitset();
  GPOperationVector new_operations;
  IterateOverRootwardEdges(
      node, [this, node, &new_operations, subsplit](
                const bool rotated, const SubsplitDAGNode *parent_node) {
        const auto parent_subsplit = parent_node->GetBitset(rotated);
        const auto gpcsp_idx = GetGPCSPIndex(parent_subsplit, subsplit);

        new_operations.push_back(IncrementWithWeightedEvolvedPLV{
            GetPLVIndex(PLVType::R_HAT, node->Id()), gpcsp_idx,
            GetPLVIndex(RPLVType(rotated), parent_node->Id())});
      });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeSBNParametersForASubsplit(const Bitset &subsplit,
                                              GPOperationVector &operations) const {
  if (parent_to_range_.count(subsplit) > 0) {
    const auto param_range = parent_to_range_.at(subsplit);
    if (param_range.second - param_range.first > 1) {
      operations.push_back(
          UpdateSBNProbabilities{param_range.first, param_range.second});
    }
  }
}

void GPDAG::UpdateRHat(size_t node_id, bool rotated,
                       GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  PLVType src_plv_type = rotated ? PLVType::R_TILDE : PLVType::R;
  const auto parent_nodes =
      rotated ? node->GetRootwardRotated() : node->GetRootwardSorted();
  GPOperationVector new_operations;
  for (size_t parent_id : parent_nodes) {
    const auto parent_node = GetDagNode(parent_id);
    auto parent_subsplit =
        rotated ? parent_node->GetBitset().RotateSubsplit() : parent_node->GetBitset();
    size_t gpcsp_idx = GetGPCSPIndex(parent_subsplit, node->GetBitset());
    new_operations.push_back(
        IncrementWithWeightedEvolvedPLV{GetPLVIndex(PLVType::R_HAT, node_id), gpcsp_idx,
                                        GetPLVIndex(src_plv_type, parent_id)});
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::UpdatePHatComputeLikelihood(size_t node_id, size_t child_node_id,
                                        bool rotated,
                                        GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  const auto child_node = GetDagNode(child_node_id);
  auto parent_subsplit =
      rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  size_t gpcsp_idx = GetGPCSPIndex(parent_subsplit, child_node->GetBitset());
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  new_operations.push_back(Likelihood{gpcsp_idx,
                                      GetPLVIndex(RPLVType(rotated), node->Id()),
                                      GetPLVIndex(PLVType::P, child_node->Id())});
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeBranchLengthUpdatePHat(size_t node_id, size_t child_node_id,
                                           bool rotated,
                                           GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  const auto child_node = GetDagNode(child_node_id);
  auto parent_subsplit =
      rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  size_t gpcsp_idx = GetGPCSPIndex(parent_subsplit, child_node->GetBitset());
  operations.push_back(OptimizeBranchLength{GetPLVIndex(PLVType::P, child_node_id),
                                            GetPLVIndex(RPLVType(rotated), node_id),
                                            gpcsp_idx});
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}
