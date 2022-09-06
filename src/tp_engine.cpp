// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "tp_engine.hpp"
#include "gp_engine.hpp"
#include "sbn_maps.hpp"

TPEngine::TPEngine(GPDAG &dag, SitePattern &site_pattern,
                   const std::string &mmap_likelihood_path, bool using_likelihoods,
                   const std::string &mmap_parsimony_path, bool using_parsimony)
    : dag_(dag),
      site_pattern_(site_pattern),
      likelihood_pvs_(mmap_likelihood_path, 0, site_pattern_.PatternCount(), 2.0),
      using_likelihoods_(using_likelihoods),
      sankoff_handler_(site_pattern, mmap_parsimony_path),
      parsimony_pvs_(sankoff_handler_.GetPSVHandler()),
      using_parsimony_(using_parsimony),
      choice_map_(dag) {
  // Initialize site pattern-based data.
  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ = EigenVectorXdOfStdVectorDouble(weights);
  // Initialize node-based data
  GrowNodeData(dag_.NodeCount(), std::nullopt, std::nullopt, true);
  // Initialize edge-based data
  GrowEdgeData(dag_.EdgeCountWithLeafSubsplits(), std::nullopt, std::nullopt, true);
  // Initialize scores.
  InitializeChoiceMap();
}

// ** General Scoring

void TPEngine::InitializeChoiceMap() { choice_map_.SelectFirstEdge(); }

Node::NodePtr TPEngine::GetTopTreeTopologyWithEdge(const EdgeId edge_id) const {
  return choice_map_.ExtractTopology(edge_id);
}

double TPEngine::GetTopTreeScore(const EdgeId edge_id) const {
  // !ISSUE #440
  Failwith("Currently no implementation.");
  return 0.0;
}

// ** Scoring by Likelihood

double TPEngine::GetTopTreeLikelihoodWithEdge(const EdgeId edge_id) {
  return top_tree_log_likelihoods_per_edge_[edge_id.value_];
}

void TPEngine::InitializeLikelihood() {
  auto &pvs = likelihood_pvs_;
  // Set all PVs to Zero
  for (NodeId node_id = 0; node_id < pvs.GetNodeCount(); node_id++) {
    for (const auto pv_type : PLVTypeEnum::Iterator()) {
      pvs.GetPV(pv_type, node_id).setZero();
    }
  }
  // Populate Leaves with Site Patterns.
  PopulateLeafLikelihoodPVsWithSitePatterns();
  // Populate Rootsplit with Stationary Distribution.
  PopulateRootLikelihoodPVsWithStationaryDistribution();
  // Rootward Pass (populate P PVs)
  const auto rootward_node_ids = dag_.RootwardNodeTraversalTrace(false);
  for (const auto node_id : rootward_node_ids) {
    PopulateRootwardLikelihoodPVForNode(node_id);
  }
  // Leafward Pass (populate R PVs)
  const auto leafward_node_ids = dag_.LeafwardNodeTraversalTrace(false);
  for (const auto node_id : leafward_node_ids) {
    PopulateLeafwardLikelihoodPVForNode(node_id);
  }
}

void TPEngine::UpdateLikelihoodsAfterDAGAddNodePair(const NNIOperation &post_nni,
                                                    const NNIOperation &pre_nni,
                                                    std::optional<size_t> new_tree_id) {
  // Copy over branch lengths.
  CopyOverEdgeDataFromPreNNIToPostNNI(post_nni, pre_nni, new_tree_id);
  // Populate PLVs.
  InitializeLikelihood();
}

void TPEngine::CopyOverEdgeDataFromPreNNIToPostNNI(const NNIOperation &post_nni,
                                                   const NNIOperation &pre_nni,
                                                   std::optional<size_t> new_tree_id) {
  new_tree_id = (new_tree_id.has_value()) ? new_tree_id.value() : GetInputTreeCount();
  input_tree_count_ = new_tree_id.value() + 1;
  // Copy over all adjacent branch lengths from those
  auto CopyBranchLengthFromCommonAdjacentNodes =
      [this, new_tree_id](const NodeId pre_node_id, const NodeId post_node_id,
                          const Direction direction, const SubsplitClade clade) {
        const auto pre_node = dag_.GetDAGNode(pre_node_id);
        for (const auto parent_id : pre_node.GetNeighbors(direction, clade)) {
          const auto pre_edge_id = dag_.GetEdgeIdx(parent_id, pre_node_id);
          const auto post_edge_id = dag_.GetEdgeIdx(parent_id, post_node_id);
          branch_lengths_[post_edge_id.value_] = branch_lengths_[pre_edge_id.value_];
          tree_source_[post_edge_id.value_] = new_tree_id.value();
        }
      };
  const auto pre_parent_id = dag_.GetDAGNodeId(pre_nni.GetParent());
  const auto pre_child_id = dag_.GetDAGNodeId(pre_nni.GetChild());
  const auto pre_edge_id = dag_.GetEdgeIdx(pre_parent_id, pre_child_id);
  const auto pre_edge = dag_.GetDAGEdge(pre_edge_id);
  const auto post_parent_id = dag_.GetDAGNodeId(post_nni.GetParent());
  const auto post_child_id = dag_.GetDAGNodeId(post_nni.GetChild());
  const auto post_edge_id = dag_.GetEdgeIdx(post_parent_id, post_child_id);
  // Copy over central edge.
  branch_lengths_[post_edge_id.value_] = branch_lengths_[pre_edge_id.value_];
  // Copy over parent and sister edges.
  CopyBranchLengthFromCommonAdjacentNodes(pre_parent_id, post_parent_id,
                                          Direction::Rootward, SubsplitClade::Left);
  CopyBranchLengthFromCommonAdjacentNodes(pre_parent_id, post_parent_id,
                                          Direction::Rootward, SubsplitClade::Right);
  CopyBranchLengthFromCommonAdjacentNodes(
      pre_parent_id, post_child_id, Direction::Leafward,
      Bitset::Opposite(pre_edge.GetSubsplitClade()));
  // Copy over left and right edges.
  NodeId post_leftchild_id;
  NodeId post_rightchild_id;
  if (pre_nni.GetSisterClade() == post_nni.GetLeftChildClade()) {
    // If post_nni swapped pre_nni sister with pre_nni left child.
    post_leftchild_id = post_parent_id;
    post_rightchild_id = post_child_id;
  } else {
    // If post_nni swapped pre_nni sister swapped with pre_nni right child.
    post_leftchild_id = post_child_id;
    post_rightchild_id = post_parent_id;
  }
  CopyBranchLengthFromCommonAdjacentNodes(pre_child_id, post_leftchild_id,
                                          Direction::Leafward, SubsplitClade::Left);
  CopyBranchLengthFromCommonAdjacentNodes(pre_child_id, post_rightchild_id,
                                          Direction::Leafward, SubsplitClade::Right);
}

double TPEngine::GetTopTreeLikelihoodWithProposedNNI(const NNIOperation &post_nni,
                                                     const NNIOperation &pre_nni) {
  using NNIClade = NNIOperation::NNIClade;
  using NNICladeEnum = NNIOperation::NNICladeEnum;
  auto &pvs = likelihood_pvs_;

  // Node ids from pre-NNI in DAG.
  NNICladeEnum::Array<PVId> pre_id_map;
  const auto pre_parent_id = dag_.GetDAGNodeId(pre_nni.GetParent());
  const auto pre_child_id = dag_.GetDAGNodeId(pre_nni.GetChild());
  const auto pre_edge_id = dag_.GetEdgeIdx(pre_parent_id, pre_child_id);
  const auto pre_sister_clade =
      Bitset::Opposite(dag_.GetDAGEdge(pre_edge_id).GetSubsplitClade());
  // PLV ids from pre-NNI in DAG.
  pre_id_map[NNIClade::ParentFocal] = pvs.GetPVIndex(PLVType::RHat, pre_parent_id);
  pre_id_map[NNIClade::ParentSister] =
      pvs.GetPVIndex(pvs.PPLVType(pre_sister_clade), pre_parent_id);
  pre_id_map[NNIClade::ChildLeft] = pvs.GetPVIndex(PLVType::PHatLeft, pre_child_id);
  pre_id_map[NNIClade::ChildRight] = pvs.GetPVIndex(PLVType::PHatRight, pre_child_id);

  // Compute P-PLV of post-child node.
  auto child_phat_pvid = likelihood_pvs_.GetSparePVIndex(0);
  // Compute R-PLV of post-parent node.
  auto parent_rfocal_pvid = likelihood_pvs_.GetSparePVIndex(1);
  // Get temp central edge.
  auto post_edge_id = EdgeId(dag_.EdgeCountWithLeafSubsplits() + 1);

  const auto clade_map =
      NNIOperation::BuildNNICladeMapFromPreNNIToNNI(pre_nni, post_nni);

  // Evolve post-leftchild with post-rightchild together.
  MultiplyPVs(child_phat_pvid, pre_id_map[clade_map[NNIClade::ChildLeft]],
              pre_id_map[clade_map[NNIClade::ChildRight]]);
  // Evolve post-parent with post-sister together.
  MultiplyPVs(parent_rfocal_pvid, pre_id_map[NNIClade::ParentFocal],
              pre_id_map[clade_map[NNIClade::ParentSister]]);

  // Get likelihood by evolving up from child to parent.
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(pre_edge_id.value_));
  PreparePerPatternLogLikelihoodsForEdge(parent_rfocal_pvid, child_phat_pvid);
  log_likelihoods_.row(post_edge_id.value_) = per_pattern_log_likelihoods_;
  auto top_tree_log_likelihood_ =
      log_likelihoods_.block(post_edge_id.value_, 0, 1, log_likelihoods_.cols()) *
      site_pattern_weights_;
  return top_tree_log_likelihood_[0];
}

// For rootward traversal. Compute the likelihood PV for a given node, by using
// the left and right child edges from the choice map of the edge below node.
void TPEngine::PopulateRootwardLikelihoodPVForNode(const NodeId node_id) {
  auto &pvs = likelihood_pvs_;
  // Iterate over all rootward edges to find the edge from the best source tree.
  EdgeId best_edge_id = FindBestEdgeAdjacentToNode(node_id, Direction::Rootward);
  const auto edge_choice = choice_map_.GetEdgeChoice(best_edge_id);
  // Evolve up from left child.
  const EdgeId left_child_edge_id = edge_choice.left_child_edge_id;
  if (left_child_edge_id != NoId) {
    EvolveLikelihoodPPVUpEdge(left_child_edge_id);
  }
  // Evolve up from right child.
  const EdgeId right_child_edge_id = edge_choice.right_child_edge_id;
  if (right_child_edge_id != NoId) {
    EvolveLikelihoodPPVUpEdge(right_child_edge_id);
  }
  // Update P-PLV from PHatLeft and PHatRight.
  const PVId pv_p = pvs.GetPVIndex(PLVType::P, node_id);
  const PVId pv_phatleft = pvs.GetPVIndex(PLVType::PHatLeft, node_id);
  const PVId pv_phatright = pvs.GetPVIndex(PLVType::PHatRight, node_id);
  if (right_child_edge_id != NoId && left_child_edge_id != NoId) {
    MultiplyPVs(pv_p, pv_phatleft, pv_phatright);
  } else if (right_child_edge_id != NoId) {
    TakePVValue(pv_p, pv_phatright);
  } else if (left_child_edge_id != NoId) {
    TakePVValue(pv_p, pv_phatleft);
  }
}

// For leafward traversal. Compute the likelihood PV for a given node, by using the
// parent and sister edges from the choice map of the edge above node.
void TPEngine::PopulateLeafwardLikelihoodPVForNode(const NodeId node_id) {
  auto &pvs = likelihood_pvs_;
  // Evolve down parent.
  // If node is not a leaf, find best edge from leafward edge's choicemap.
  EdgeId best_edge_id = EdgeId(NoId);
  if (!dag_.IsNodeLeaf(node_id)) {
    best_edge_id = FindBestEdgeAdjacentToNode(node_id, Direction::Leafward);
    const auto parent_edge_id = choice_map_.GetEdgeChoice(best_edge_id).parent_edge_id;
    if (parent_edge_id != NoId) {
      const auto parent_edge = dag_.GetDAGEdge(parent_edge_id);
      if (parent_edge.GetParent() != dag_.GetDAGRootNodeId()) {
        EvolveLikelihoodRPVDownEdge(parent_edge_id);
      }
    }
  }
  // If node is a leaf, find best edge from rootward edge's choicemap.
  else {
    best_edge_id = FindBestEdgeAdjacentToNode(node_id, Direction::Rootward);
    const auto edge = dag_.GetDAGEdge(best_edge_id);
    best_edge_id = choice_map_.GetEdgeChoice(best_edge_id).parent_edge_id;
    best_edge_id = (edge.GetSubsplitClade() == SubsplitClade::Left)
                       ? choice_map_.GetEdgeChoice(best_edge_id).left_child_edge_id
                       : choice_map_.GetEdgeChoice(best_edge_id).right_child_edge_id;
    EvolveLikelihoodRPVDownEdge(best_edge_id);
  }

  // Evolve with sister.
  const PVId pv_center_index = pvs.GetPVIndex(PLVType::RHat, node_id);
  const PVId pv_left_index = pvs.GetPVIndex(PLVType::RLeft, node_id);
  const PVId pv_right_index = pvs.GetPVIndex(PLVType::RRight, node_id);
  const PVId pv_left_child_index = pvs.GetPVIndex(PLVType::PHatLeft, node_id);
  const PVId pv_right_child_index = pvs.GetPVIndex(PLVType::PHatRight, node_id);
  MultiplyPVs(pv_left_index, pv_right_child_index, pv_center_index);
  MultiplyPVs(pv_right_index, pv_left_child_index, pv_center_index);
}

void TPEngine::PopulateLeafLikelihoodPVsWithSitePatterns() {
  auto &pvs = likelihood_pvs_;
  for (PVId pv_id = 0; pv_id < pvs.GetPVCount(); pv_id++) {
    pvs.GetPV(pv_id).setZero();
  }
  NodeId taxon_idx = 0;
  for (const auto &pattern : site_pattern_.GetPatterns()) {
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative symbol!");
      for (const auto pv_type : {PLVType::P}) {
        if (symbol == MmappedNucleotidePLV::base_count_) {  // Gap character.
          pvs.GetPV(pvs.GetPVIndex(pv_type, taxon_idx)).col(site_idx).setConstant(1.);
        } else if (symbol < MmappedNucleotidePLV::base_count_) {
          pvs.GetPV(pvs.GetPVIndex(pv_type, taxon_idx))(symbol, site_idx) = 1.;
        }
      }
      site_idx++;
    }
    taxon_idx++;
  }
}

void TPEngine::PopulateRootLikelihoodPVsWithStationaryDistribution() {
  for (const auto pv_type : {PLVType::RHat, PLVType::RLeft, PLVType::RRight}) {
    auto &pv = likelihood_pvs_.GetPV(pv_type, dag_.GetDAGRootNodeId());
    for (Eigen::Index row_idx = 0; row_idx < pv.rows(); ++row_idx) {
      pv.row(row_idx).array() = stationary_distribution_(row_idx);
    }
  }
  for (const auto pv_type : {PLVType::RHat}) {
    for (const auto &node_id : dag_.GetRootsplitNodeIds()) {
      auto &pv = likelihood_pvs_.GetPV(pv_type, node_id);
      for (Eigen::Index row_idx = 0; row_idx < pv.rows(); ++row_idx) {
        pv.row(row_idx).array() = stationary_distribution_(row_idx);
      }
    }
  }
}

void TPEngine::ComputeLikelihoods() {
  for (NodeId node_id = 0; node_id < dag_.NodeCount(); node_id++) {
    const auto node = dag_.GetDAGNode(node_id);
    dag_.IterateOverLeafwardEdges(
        node, [this, node](const bool is_edge_on_left, SubsplitDAGNode child_node) {
          auto &pvs = likelihood_pvs_;
          const EdgeId edge_id = dag_.GetEdgeIdx(node.Id(), child_node.Id());
          const PVId parent_pvid =
              pvs.GetPVIndex(PLVHandler::RPLVType(is_edge_on_left), node.Id());
          const PVId child_pvid = pvs.GetPVIndex(PLVType::P, child_node.Id());
          ComputeLikelihood(edge_id, parent_pvid, child_pvid);
        });
  };

  top_tree_log_likelihoods_per_edge_ =
      log_likelihoods_.block(0, 0, dag_.EdgeCountWithLeafSubsplits(),
                             log_likelihoods_.cols()) *
      site_pattern_weights_;
}

EdgeId TPEngine::FindBestEdgeAdjacentToNode(const NodeId node_id,
                                            const Direction direction) const {
  const auto node = dag_.GetDAGNode(node_id);
  EdgeId best_edge_id = EdgeId(NoId);
  size_t best_tree_source = GetInputTreeCount();
  for (const auto focal_clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto rootward_node_id : node.GetNeighbors(direction, focal_clade)) {
      const auto edge_id = dag_.GetEdgeIdx(rootward_node_id, node_id);
      if (best_tree_source > tree_source_[edge_id.value_]) {
        best_tree_source = tree_source_[edge_id.value_];
        best_edge_id = edge_id;
      }
    }
  }
  return best_edge_id;
}

EdgeId TPEngine::FindBestEdgeAdjacentToNode(const NodeId node_id,
                                            const Direction direction,
                                            const SubsplitClade clade) const {
  const auto node = dag_.GetDAGNode(node_id);
  EdgeId best_edge_id = EdgeId(NoId);
  size_t best_tree_source = GetInputTreeCount();
  for (const auto rootward_node_id : node.GetNeighbors(direction, clade)) {
    const auto edge_id = dag_.GetEdgeIdx(rootward_node_id, node_id);
    if (best_tree_source > tree_source_[edge_id.value_]) {
      best_tree_source = tree_source_[edge_id.value_];
      best_edge_id = edge_id;
    }
  }
  return best_edge_id;
}

void TPEngine::EvolveLikelihoodPPVUpEdge(const EdgeId edge_id) {
  auto &pvs = likelihood_pvs_;
  const auto edge = dag_.GetDAGEdge(edge_id);
  const PVId pv_parent_index =
      pvs.GetPVIndex(PLVHandler::PPLVType(edge.GetSubsplitClade()), edge.GetParent());
  const PVId pv_child_index = pvs.GetPVIndex(PLVType::P, edge.GetChild());
  SetToEvolvedPV(pv_parent_index, edge_id, pv_child_index);
}

void TPEngine::EvolveLikelihoodRPVDownEdge(const EdgeId edge_id) {
  auto &pvs = likelihood_pvs_;
  const auto edge = dag_.GetDAGEdge(edge_id);
  const PLVType parent_focal_rpv = PLVHandler::RPLVType(edge.GetSubsplitClade());
  const PVId pv_parent_index = pvs.GetPVIndex(parent_focal_rpv, edge.GetParent());
  const PVId pv_child_index = pvs.GetPVIndex(PLVType::RHat, edge.GetChild());
  SetToEvolvedPV(pv_child_index, edge_id, pv_parent_index);
}

// ** Scoring by Parsimony

double TPEngine::GetTopTreeParsimonyWithEdge(const EdgeId edge_id) {
  return top_tree_parsimony_per_edge_[edge_id.value_];
}

void TPEngine::PopulateLeafParsimonyPVsWithSitePatterns() {
  auto &pvs = parsimony_pvs_;
  for (auto &pv : pvs.GetPVs()) {
    pv.setZero();
  }
  NodeId taxon_idx = 0;
  for (const auto &pattern : site_pattern_.GetPatterns()) {
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative symbol!");
      for (const auto pv_type : {PSVType::PLeft, PSVType::PRight}) {
        if (symbol == MmappedNucleotidePLV::base_count_) {  // Gap character.
          pvs.GetPV(pvs.GetPVIndex(pv_type, taxon_idx)).col(site_idx).setConstant(1.);
        } else if (symbol < MmappedNucleotidePLV::base_count_) {
          pvs.GetPV(pvs.GetPVIndex(pv_type, taxon_idx))(symbol, site_idx) = 1.;
        }
      }
      site_idx++;
    }
    taxon_idx++;
  }
}

void TPEngine::InitializeParsimony() {
  auto &pvs = parsimony_pvs_;
  sankoff_handler_.Resize(dag_.NodeCount());
  // Set all PVs to Zero
  for (NodeId node_id = 0; node_id < pvs.GetNodeCount(); node_id++) {
    for (const auto pv_type : PSVTypeEnum::Iterator()) {
      pvs.GetPV(pv_type, node_id).setZero();
    }
  }
  // Populate Leaves with Site Patterns.
  sankoff_handler_.GenerateLeafPartials();

  // Rootward Pass (populate P PVs)
  const auto rootward_node_ids = dag_.RootwardNodeTraversalTrace(false);
  for (const auto node_id : rootward_node_ids) {
    PopulateRootwardParsimonyPVForNode(node_id);
  }
  // Leafward Pass (populate R PVs)
  const auto leafward_node_ids = dag_.LeafwardNodeTraversalTrace(false);
  for (const auto node_id : leafward_node_ids) {
    PopulateLeafwardParsimonyPVForNode(node_id);
  }
}

void TPEngine::PopulateRootwardParsimonyPVForNode(const NodeId node_id) {
  // Iterate over all rootward edges to find the edge from the best source tree.
  EdgeId best_edge_id = FindBestEdgeAdjacentToNode(node_id, Direction::Rootward);
  if (!dag_.IsNodeLeaf(node_id)) {
    const auto edge_choice = choice_map_.GetEdgeChoice(best_edge_id);
    // Accumulate parsimony from left and right child.
    const EdgeId left_child_edge_id = edge_choice.left_child_edge_id;
    const NodeId left_child_node_id = dag_.GetDAGEdge(left_child_edge_id).GetChild();
    const EdgeId right_child_edge_id = edge_choice.right_child_edge_id;
    const NodeId right_child_node_id = dag_.GetDAGEdge(right_child_edge_id).GetChild();

    sankoff_handler_.PopulateRootwardParsimonyPVForNode(node_id, left_child_node_id,
                                                        right_child_node_id);
  }
}

void TPEngine::PopulateLeafwardParsimonyPVForNode(const NodeId node_id) {
  // Evolve down parent.
  // If node is not a leaf, find best edge from leafward edge's choicemap.
  EdgeId best_edge_id = EdgeId(NoId);
  if (!dag_.IsNodeLeaf(node_id)) {
    best_edge_id = FindBestEdgeAdjacentToNode(node_id, Direction::Leafward);
    const auto best_edge = dag_.GetDAGEdge(best_edge_id);
    const NodeId focal_node_id = best_edge.GetChild();
    const EdgeId sister_edge_id =
        choice_map_.GetEdgeChoice(best_edge_id).sister_edge_id;
    const NodeId sister_node_id = dag_.GetDAGEdge(sister_edge_id).GetChild();
    const NodeId left_child_node_id =
        (best_edge.GetSubsplitClade() == SubsplitClade::Left) ? focal_node_id
                                                              : sister_node_id;
    const NodeId right_child_node_id =
        (best_edge.GetSubsplitClade() == SubsplitClade::Left) ? sister_node_id
                                                              : focal_node_id;

    sankoff_handler_.PopulateLeafwardParsimonyPVForNode(node_id, left_child_node_id,
                                                        right_child_node_id);
  }
}

void TPEngine::ComputeParsimonies() {
  for (EdgeId edge_id = 0; edge_id < dag_.EdgeCountWithLeafSubsplits(); edge_id++) {
    NodeId node_id = dag_.GetDAGEdge(edge_id).GetChild();
    top_tree_parsimony_per_edge_[node_id.value_] =
        sankoff_handler_.ParsimonyScore(node_id);
  }
}

void TPEngine::UpdateDAGAfterAddNodePairByParsimony(const NNIOperation &nni_op) {
  // !ISSUE #440
  Failwith("Currently no implementation.");
}

// ** Partial Vector Operations

void TPEngine::TakePVValue(const PVId dest_id, const PVId src_id) {
  auto &pvs = likelihood_pvs_;
  pvs.GetPV(dest_id).array() = pvs.GetPV(src_id).array();
}

void TPEngine::MultiplyPVs(const PVId dest_id, const PVId src1_id, const PVId src2_id) {
  auto &pvs = likelihood_pvs_;
  pvs.GetPV(dest_id).array() = pvs.GetPV(src1_id).array() * pvs.GetPV(src2_id).array();
  // #462: Need to add rescaling to PVs.
}

void TPEngine::ComputeLikelihood(const EdgeId dest_id, const PVId child_id,
                                 const PVId parent_id) {
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(dest_id.value_));
  PreparePerPatternLogLikelihoodsForEdge(parent_id, child_id);
  log_likelihoods_.row(dest_id.value_) = per_pattern_log_likelihoods_;
}

void TPEngine::SetToEvolvedPV(const PVId dest_id, const EdgeId edge_id,
                              const PVId src_id) {
  auto &pvs = likelihood_pvs_;
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(edge_id.value_));
  pvs.GetPV(dest_id).array() = (transition_matrix_ * pvs.GetPV(src_id));
}

void TPEngine::MultiplyWithEvolvedPV(const PVId dest_id, const EdgeId edge_id,
                                     const PVId src_id) {
  auto &pvs = likelihood_pvs_;
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(edge_id.value_));
  pvs.GetPV(dest_id).array() =
      pvs.GetPV(dest_id).array() * (transition_matrix_ * pvs.GetPV(src_id)).array();
}

void TPEngine::SetTransitionMatrixToHaveBranchLength(const double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void TPEngine::PreparePerPatternLogLikelihoodsForEdge(const PVId src1_id,
                                                      const PVId src2_id) {
  auto &pvs = likelihood_pvs_;
  per_pattern_log_likelihoods_ =
      (pvs.GetPV(src1_id).transpose() * transition_matrix_ * pvs.GetPV(src2_id))
          .diagonal()
          .array()
          .log();
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
    log_likelihoods_.conservativeResize(GetAllocatedEdgeCount(),
                                        site_pattern_.PatternCount());
    top_tree_log_likelihoods_per_edge_.conservativeResize(GetAllocatedEdgeCount());
    top_tree_parsimony_per_edge_.conservativeResize(GetAllocatedEdgeCount());
  }
  // Resize to fit without deallocating unused memory.
  branch_lengths_.conservativeResize(GetPaddedEdgeCount());
  log_likelihoods_.conservativeResize(GetPaddedEdgeCount(),
                                      site_pattern_.PatternCount());
  top_tree_log_likelihoods_per_edge_.conservativeResize(GetPaddedEdgeCount());
  top_tree_parsimony_per_edge_.conservativeResize(GetPaddedEdgeCount());

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
  Reindexer::ReindexInPlace<EigenVectorXd, double>(top_tree_log_likelihoods_per_edge_,
                                                   edge_reindexer, GetEdgeCount());
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

// ** Tree Collection

void TPEngine::SetChoiceMapByTakingFirst(const RootedTreeCollection &tree_collection,
                                         const BitsetSizeMap &edge_indexer) {
  input_tree_count_ = tree_collection.TreeCount();
  const size_t tree_id_max = GetInputTreeCount();
  tree_source_.resize(GetEdgeCount(), tree_id_max);
  std::fill(tree_source_.begin(), tree_source_.end(), tree_id_max);
  // Set tree source map for each edge in DAG.
  auto set_tree_source = [tree_id_max, this](
                             const EdgeId edge_idx, const Bitset &edge_bitset,
                             const RootedTree &tree, const size_t tree_id,
                             const Node *focal_node) {
    if (tree_source_[edge_idx.value_] == tree_id_max) {
      tree_source_[edge_idx.value_] = tree_id;
    }
  };
  RootedSBNMaps::FunctionOverRootedTreeCollection(set_tree_source, tree_collection,
                                                  edge_indexer, branch_lengths_.size());
  // Set tree source map for rootsplit edges from the best tree.
  const auto root_node = dag_.GetDAGNode(dag_.GetDAGRootNodeId());
  for (const auto rootsplit_node_id : dag_.GetRootsplitNodeIds()) {
    const auto rootsplit_node = dag_.GetDAGNode(rootsplit_node_id);
    const auto rootsplit_edge_id = dag_.GetEdgeIdx(root_node.Id(), rootsplit_node.Id());
    size_t best_tree_source = tree_id_max;
    dag_.IterateOverLeafwardEdges(
        rootsplit_node, [this, &rootsplit_node, &rootsplit_edge_id, &best_tree_source](
                            bool is_edge_on_left, SubsplitDAGNode child_node) {
          const auto edge_id = dag_.GetEdgeIdx(rootsplit_node.Id(), child_node.Id());
          if (best_tree_source > tree_source_[edge_id.value_]) {
            best_tree_source = tree_source_[edge_id.value_];
            tree_source_[rootsplit_edge_id.value_] = best_tree_source;
          }
        });
  }
  // Build choice map from tree source.
  for (EdgeId edge_id = 0; edge_id < GetEdgeCount(); edge_id++) {
    const auto edge = dag_.GetDAGEdge(edge_id);
    EdgeId best_edge_id;
    // Get parent EdgeChoice.
    best_edge_id = FindBestEdgeAdjacentToNode(edge.GetParent(), Direction::Rootward);
    choice_map_.SetEdgeChoice(edge_id, ChoiceMap::AdjacentEdge::Parent, best_edge_id);
    // Get sister EdgeChoice.
    best_edge_id =
        FindBestEdgeAdjacentToNode(edge.GetParent(), Direction::Leafward,
                                   Bitset::Opposite(edge.GetSubsplitClade()));
    choice_map_.SetEdgeChoice(edge_id, ChoiceMap::AdjacentEdge::Sister, best_edge_id);
    // Get left child EdgeChoice.
    best_edge_id = FindBestEdgeAdjacentToNode(edge.GetChild(), Direction::Leafward,
                                              SubsplitClade::Left);
    choice_map_.SetEdgeChoice(edge_id, ChoiceMap::AdjacentEdge::LeftChild,
                              best_edge_id);
    // Get right child EdgeChoice.
    best_edge_id = FindBestEdgeAdjacentToNode(edge.GetChild(), Direction::Leafward,
                                              SubsplitClade::Right);
    choice_map_.SetEdgeChoice(edge_id, ChoiceMap::AdjacentEdge::RightChild,
                              best_edge_id);
  }
}

void TPEngine::SetBranchLengthByTakingFirst(const RootedTreeCollection &tree_collection,
                                            const BitsetSizeMap &edge_indexer) {
  // Unique edges in collection should be the same as the number of total edges in DAG
  // created from collection.
  branch_lengths_.setZero();
  EigenVectorXi observed_edge_counts = EigenVectorXi::Zero(GetEdgeCount());
  // Set branch lengths on first occurance.
  auto set_first_branch_length = [&observed_edge_counts, this](
                                     const EdgeId edge_idx, const Bitset &edge_bitset,
                                     const RootedTree &tree, const size_t tree_id,
                                     const Node *focal_node) {
    if (observed_edge_counts(edge_idx.value_) == 0) {
      branch_lengths_(edge_idx.value_) = tree.BranchLength(focal_node);
      observed_edge_counts(edge_idx.value_)++;
    }
  };
  RootedSBNMaps::FunctionOverRootedTreeCollection(
      set_first_branch_length, tree_collection, edge_indexer, branch_lengths_.size());
  // Either set branch length to the first occurance if edge exists in collection,
  // otherwise set length to default.
  for (EdgeId edge_idx = 0; edge_idx < GetEdgeCount(); edge_idx++) {
    if (observed_edge_counts(edge_idx.value_) == 0) {
      branch_lengths_(edge_idx.value_) = default_branch_length_;
    }
  }
}

// ** I/O

std::string TPEngine::LikelihoodPVToString(const PVId pv_id) const {
  return likelihood_pvs_.ToString(pv_id);
}

std::string TPEngine::LogLikelihoodMatrixToString() const {
  std::stringstream out;
  for (Eigen::Index i = 0; i < log_likelihoods_.rows(); i++) {
    for (Eigen::Index j = 0; j < log_likelihoods_.cols(); j++) {
      out << "[" << i << "," << j << "]: " << log_likelihoods_(i, j) << "\t";
    }
    out << std::endl;
  }
  return out.str();
}

std::string TPEngine::ParsimonyPVToString(const PVId pv_id) const {
  return parsimony_pvs_.ToString(pv_id);
}
