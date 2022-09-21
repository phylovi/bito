// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "gp_dag.hpp"

#include "numerical_utils.hpp"

using namespace GPOperations;  // NOLINT
using PLVType = PLVNodeHandler::PLVType;

size_t GPDAG::GetPLVIndex(PLVType plv_type, NodeId node_id) const {
  return PLVNodeHandler::GetPVIndex(plv_type, node_id, NodeCountWithoutDAGRoot())
      .value_;
}

// The R PLV update that corresponds to our rotation status.
GPOperation GPDAG::RUpdateOfRotated(NodeId node_id, bool is_edge_on_left) const {
  return is_edge_on_left ? Multiply{GetPLVIndex(PLVType::RLeft, node_id),
                                    GetPLVIndex(PLVType::RHat, node_id),
                                    GetPLVIndex(PLVType::PHatRight, node_id)}
                         : Multiply{GetPLVIndex(PLVType::RRight, node_id),
                                    GetPLVIndex(PLVType::RHat, node_id),
                                    GetPLVIndex(PLVType::PHatLeft, node_id)};
}

// After this traversal, we will have optimized branch lengths, but we cannot assume
// that all of the PLVs are in a valid state.
//
// Update the terminology in this function as part of #288.
GPOperationVector GPDAG::ApproximateBranchLengthOptimization() const {
  GPOperationVector operations;
  SubsplitDAG::DepthFirstWithAction(
      GetRootsplitNodeIds(),
      SubsplitDAGTraversalAction(
          // BeforeNode
          [this, &operations](NodeId node_id) {
            if (!GetDAGNode(node_id).IsRootsplit()) {
              // Update R-hat if we're not at the root.
              UpdateRHat(node_id, operations);
            }
          },
          // AfterNode
          [this, &operations](NodeId node_id) {
            // Make P the elementwise product ("o") of the two P PLVs for the
            // node-clades.
            operations.push_back(Multiply{GetPLVIndex(PLVType::P, node_id),
                                          GetPLVIndex(PLVType::PHatRight, node_id),
                                          GetPLVIndex(PLVType::PHatLeft, node_id)});
          },
          // BeforeNodeClade
          [this, &operations](NodeId node_id, bool is_edge_on_left) {
            const PLVType p_hat_plv_type =
                is_edge_on_left ? PLVType::PHatLeft : PLVType::PHatRight;
            // Update the R PLV corresponding to our rotation status.
            operations.push_back(RUpdateOfRotated(node_id, is_edge_on_left));
            // Zero out the node-clade PLV so we can fill it as part of VisitEdge.
            operations.push_back(ZeroPLV{GetPLVIndex(p_hat_plv_type, node_id)});
          },
          // VisitEdge
          [this, &operations](NodeId node_id, NodeId child_id, bool is_edge_on_left) {
            // #310 this is temporary:
            // We do a full PLV population and then marginal likelihood calculation.
            // GPOperations::AppendGPOperations(operations, PopulatePLVs());
            // GPOperations::AppendGPOperations(operations, MarginalLikelihood());

            // Optimize each branch for a given node-clade and accumulate the resulting
            // P-hat PLVs in the parent node.
            OptimizeBranchLengthUpdatePHat(node_id, child_id, is_edge_on_left,
                                           operations);
          }));
  return operations;
}

// After this traversal, we will have optimized branch lengths, but we cannot assume
// that all of the PLVs are in a valid state.
//
// Update the terminology in this function as part of #288.
GPOperationVector GPDAG::BranchLengthOptimization() {
  GPOperationVector operations;
  DepthFirstWithTidyAction(
      GetRootsplitNodeIds(),
      TidySubsplitDAGTraversalAction(
          // BeforeNode
          [this, &operations](NodeId node_id) {
            if (!GetDAGNode(node_id).IsRootsplit()) {
              // Update R-hat if we're not at the root.
              UpdateRHat(node_id, operations);
            }
          },
          // AfterNode
          [this, &operations](NodeId node_id) {
            // Make P the elementwise product ("o") of the two P PLVs for the
            // node-clades.
            operations.push_back(Multiply{GetPLVIndex(PLVType::P, node_id),
                                          GetPLVIndex(PLVType::PHatRight, node_id),
                                          GetPLVIndex(PLVType::PHatLeft, node_id)});
          },
          // BeforeNodeClade
          [this, &operations](NodeId node_id, bool is_edge_on_left) {
            const PLVType p_hat_plv_type =
                is_edge_on_left ? PLVType::PHatLeft : PLVType::PHatRight;
            // Update the R PLV corresponding to our rotation status.
            operations.push_back(RUpdateOfRotated(node_id, is_edge_on_left));
            // Zero out the node-clade PLV so we can fill it as part of VisitEdge.
            operations.push_back(ZeroPLV{GetPLVIndex(p_hat_plv_type, node_id)});
          },
          // ModifyEdge
          [this, &operations](NodeId node_id, NodeId child_id, bool is_edge_on_left) {
            // Optimize each branch for a given node-clade and accumulate the resulting
            // P-hat PLVs in the parent node.
            OptimizeBranchLengthUpdatePHat(node_id, child_id, is_edge_on_left,
                                           operations);
          },
          // UpdateEdge
          [this, &operations](NodeId node_id, NodeId child_id, bool is_edge_on_left) {
            // Accumulate all P-hat PLVs in the parent node without optimization.
            // #321 I don't think we need this Likelihood call... just the update PHat.
            UpdatePHatComputeLikelihood(node_id, child_id, is_edge_on_left, operations);
          }));
  return operations;
}

GPOperationVector GPDAG::ComputeLikelihoods() const {
  GPOperationVector operations;
  IterateOverRealNodes([this, &operations](SubsplitDAGNode node) {
    IterateOverLeafwardEdges(
        node, [this, node, &operations](const bool is_edge_on_left,
                                        SubsplitDAGNode child_node) {
          const auto gpcsp_idx = GetEdgeIdx(node.Id(), child_node.Id());
          operations.push_back(Likelihood{
              gpcsp_idx.value_,
              GetPLVIndex(PLVNodeHandler::RPLVType(is_edge_on_left), node.Id()),
              GetPLVIndex(PLVType::P, child_node.Id())});
        });
  });

  const auto marginal_likelihood_operations = MarginalLikelihood();
  operations.insert(operations.end(), marginal_likelihood_operations.begin(),
                    marginal_likelihood_operations.end());

  return operations;
}

GPOperationVector GPDAG::LeafwardPass() const {
  return LeafwardPass(LeafwardNodeTraversalTrace(false));
}

GPOperationVector GPDAG::MarginalLikelihood() const {
  GPOperationVector operations = {GPOperations::ResetMarginalLikelihood{}};
  for (const auto &rootsplit_id : GetRootsplitNodeIds()) {
    operations.push_back(GPOperations::IncrementMarginalLikelihood{
        GetPLVIndex(PLVType::RHat, NodeId(rootsplit_id)),
        GetEdgeIdx(GetDAGRootNodeId(), NodeId(rootsplit_id)).value_,
        GetPLVIndex(PLVType::P, NodeId(rootsplit_id))});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass() const {
  return RootwardPass(RootwardNodeTraversalTrace(false));
}

GPOperationVector GPDAG::OptimizeSBNParameters() const {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (const auto &id : LeafwardNodeTraversalTrace(false)) {
    const auto node = GetDAGNode(id);
    OptimizeSBNParametersForASubsplit(node.GetBitset(), operations);
    OptimizeSBNParametersForASubsplit(node.GetBitset().SubsplitRotate(), operations);
  }
  operations.push_back(UpdateSBNProbabilities{0, RootsplitCount()});
  return operations;
}

GPOperationVector GPDAG::SetLeafwardZero() const {
  GPOperationVector operations;
  for (NodeId node_id = NodeId(0); node_id < NodeCountWithoutDAGRoot(); node_id++) {
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::RHat, node_id)});
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::RRight, node_id)});
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::RLeft, node_id)});
  }
  return operations;
}

GPOperationVector GPDAG::SetRhatToStationary() const {
  GPOperationVector operations;
  for (const auto &rootsplit_id : GetRootsplitNodeIds()) {
    auto rootsplit_node_id = NodeId(rootsplit_id);
    auto root_gpcsp_idx = GetEdgeIdx(GetDAGRootNodeId(), rootsplit_node_id);
    operations.push_back(SetToStationaryDistribution{
        GetPLVIndex(PLVType::RHat, rootsplit_node_id), root_gpcsp_idx.value_});
  }
  return operations;
}

GPOperationVector GPDAG::SetRootwardZero() const {
  GPOperationVector operations;
  for (NodeId node_id = NodeId(TaxonCount()); node_id < NodeCountWithoutDAGRoot();
       node_id++) {
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::P, node_id)});
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::PHatRight, node_id)});
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::PHatLeft, node_id)});
  }
  return operations;
}

GPOperationVector GPDAG::LeafwardPass(const NodeIdVector &visit_order) const {
  GPOperationVector operations;
  for (const auto node_id : visit_order) {
    auto node = GetDAGNode(node_id);
    // Build rhat(s) via rhat(s) += \sum_t q(s|t) P'(s|t) r(t)
    AddRhatOperations(node, operations);
    // Multiply to get r(s_right) = rhat(s) \circ phat(s_left).
    operations.push_back(Multiply{GetPLVIndex(PLVType::RRight, node_id),
                                  GetPLVIndex(PLVType::RHat, node_id),
                                  GetPLVIndex(PLVType::PHatLeft, node_id)});
    // Multiply to get r(s_left) = rhat(s) \circ phat(s_right).
    operations.push_back(Multiply{GetPLVIndex(PLVType::RLeft, node_id),
                                  GetPLVIndex(PLVType::RHat, node_id),
                                  GetPLVIndex(PLVType::PHatRight, node_id)});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass(const NodeIdVector &visit_order) const {
  GPOperationVector operations;
  for (const auto node_id : visit_order) {
    const auto node = GetDAGNode(node_id);
    if (!node.IsLeaf()) {
      // Build phat(s_right).
      AddPhatOperations(node, false, operations);
      // Build phat(s_left).
      AddPhatOperations(node, true, operations);
      // Multiply to get p(s) = phat(s_left) \circ phat(s_right).
      operations.push_back(Multiply{node_id.value_,
                                    GetPLVIndex(PLVType::PHatRight, node_id),
                                    GetPLVIndex(PLVType::PHatLeft, node_id)});
    }
  }
  return operations;
}

GPOperationVector GPDAG::PopulatePLVs() const {
  GPOperationVector operations;
  GPOperations::AppendGPOperations(operations, SetRootwardZero());
  GPOperations::AppendGPOperations(operations, SetLeafwardZero());
  GPOperations::AppendGPOperations(operations, SetRhatToStationary());
  GPOperations::AppendGPOperations(operations, RootwardPass());
  GPOperations::AppendGPOperations(operations, LeafwardPass());
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

void GPDAG::AddPhatOperations(SubsplitDAGNode node, bool is_edge_on_left,
                              GPOperationVector &operations) const {
  PLVType plv_type = PLVNodeHandler::PPLVType(is_edge_on_left);
  const auto parent_id = node.Id();
  const auto dest_idx = GetPLVIndex(plv_type, node.Id());
  GPOperationVector new_operations;
  for (const auto &child_id : node.GetLeafward(is_edge_on_left)) {
    const auto gpcsp_idx = GetEdgeIdx(parent_id, NodeId(child_id));
    new_operations.push_back(IncrementWithWeightedEvolvedPLV{
        dest_idx, gpcsp_idx.value_, GetPLVIndex(PLVType::P, NodeId(child_id))});
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::AddRhatOperations(SubsplitDAGNode node,
                              GPOperationVector &operations) const {
  GPOperationVector new_operations;
  IterateOverRootwardEdges(
      node, [this, node, &new_operations](const bool is_edge_on_left,
                                          SubsplitDAGNode parent_node) {
        new_operations.push_back(IncrementWithWeightedEvolvedPLV{
            GetPLVIndex(PLVType::RHat, node.Id()),
            GetEdgeIdx(parent_node.Id(), node.Id()).value_,
            GetPLVIndex(PLVNodeHandler::RPLVType(is_edge_on_left), parent_node.Id())});
      });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeSBNParametersForASubsplit(const Bitset &subsplit,
                                              GPOperationVector &operations) const {
  if (parent_to_child_range_.count(subsplit) > 0) {
    const auto &[edge_begin, edge_end] = parent_to_child_range_.at(subsplit);
    if (edge_begin.value_ - edge_end.value_ > 1) {
      operations.push_back(UpdateSBNProbabilities{edge_begin.value_, edge_end.value_});
    }
  }
}

void GPDAG::UpdateRHat(NodeId node_id, GPOperationVector &operations) const {
  operations.push_back(ZeroPLV{GetPLVIndex(PLVType::RHat, node_id)});
  GPOperationVector new_operations;
  const auto node = GetDAGNode(node_id);
  for (const bool is_edge_on_left : {false, true}) {
    PLVType src_plv_type = is_edge_on_left ? PLVType::RLeft : PLVType::RRight;
    for (auto parent_id : node.GetRootward(is_edge_on_left)) {
      new_operations.push_back(IncrementWithWeightedEvolvedPLV{
          GetPLVIndex(PLVType::RHat, node_id), GetEdgeIdx(parent_id, node_id).value_,
          GetPLVIndex(src_plv_type, parent_id)});
    }
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

// #311 there's some work to be done here.
// There's a lot of common code between this function and the next.
// Also, the prep for marginalization isn't actually working correctly: we need to
// gather more operations first.
void GPDAG::UpdatePHatComputeLikelihood(NodeId node_id, NodeId child_node_id,
                                        bool is_edge_on_left,
                                        GPOperationVector &operations) const {
  const auto gpcsp_idx = GetEdgeIdx(node_id, child_node_id);
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(is_edge_on_left ? PLVType::PHatLeft : PLVType::PHatRight, node_id),
      gpcsp_idx.value_,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  new_operations.push_back(Likelihood{
      gpcsp_idx.value_, GetPLVIndex(PLVNodeHandler::RPLVType(is_edge_on_left), node_id),
      GetPLVIndex(PLVType::P, child_node_id)});
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeBranchLengthUpdatePHat(NodeId node_id, NodeId child_node_id,
                                           bool is_edge_on_left,
                                           GPOperationVector &operations) const {
  EdgeId gpcsp_idx = GetEdgeIdx(node_id, child_node_id);
  operations.push_back(OptimizeBranchLength{
      GetPLVIndex(PLVType::P, child_node_id),
      GetPLVIndex(PLVNodeHandler::RPLVType(is_edge_on_left), node_id),
      gpcsp_idx.value_});
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(is_edge_on_left ? PLVType::PHatLeft : PLVType::PHatRight, node_id),
      gpcsp_idx.value_,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

QuartetHybridRequest GPDAG::QuartetHybridRequestOf(NodeId parent_id,
                                                   bool is_focal_on_left,
                                                   NodeId child_id) const {
  QuartetTipVector rootward_tips;
  IterateOverRootwardEdgesAndParents(
      GetDAGNode(parent_id),
      [this, &rootward_tips](const EdgeId gpcsp_idx, const bool is_rootward_on_left,
                             const NodeId grandparent_id) {
        rootward_tips.emplace_back(
            grandparent_id.value_,
            GetPLVIndex(PLVNodeHandler::RPLVType(is_rootward_on_left), grandparent_id),
            gpcsp_idx.value_);
      });

  QuartetTipVector sister_tips;
  const auto &parent_node = GetDAGNode(parent_id);
  const bool is_sister_edge_on_left = !is_focal_on_left;
  IterateOverLeafwardEdges(
      parent_node, is_sister_edge_on_left,
      [this, &parent_node, &sister_tips](SubsplitDAGNode sister_node) {
        const auto sister_id = sister_node.Id();
        sister_tips.emplace_back(
            sister_id.value_, GetPLVIndex(PLVType::P, sister_id),
            GetEdgeIdx(parent_node.GetBitset(), sister_node.GetBitset()).value_);
      });

  QuartetTipVector rotated_tips;
  QuartetTipVector sorted_tips;
  IterateOverLeafwardEdgesAndChildren(
      GetDAGNode(child_id), [this, &rotated_tips, &sorted_tips](
                                const EdgeId gpcsp_idx, const bool is_leafward_on_left,
                                const NodeId grandchild_id) {
        if (is_leafward_on_left) {
          rotated_tips.emplace_back(grandchild_id.value_,
                                    GetPLVIndex(PLVType::P, grandchild_id),
                                    gpcsp_idx.value_);
        } else {
          sorted_tips.emplace_back(grandchild_id.value_,
                                   GetPLVIndex(PLVType::P, grandchild_id),
                                   gpcsp_idx.value_);
        }
      });
  return QuartetHybridRequest(GetEdgeIdx(parent_id, child_id).value_,
                              std::move(rootward_tips), std::move(sister_tips),
                              std::move(rotated_tips), std::move(sorted_tips));
}
