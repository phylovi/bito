// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "tp_engine.hpp"
#include "gp_engine.hpp"
#include "sbn_maps.hpp"
#include "optimization.hpp"

TPEngine::TPEngine(GPDAG &dag, SitePattern &site_pattern,
                   const std::string &mmap_likelihood_path, bool using_likelihoods,
                   const std::string &mmap_parsimony_path, bool using_parsimony,
                   bool use_gradients)
    : dag_(dag),
      site_pattern_(site_pattern),
      likelihood_pvs_(mmap_likelihood_path, 0, site_pattern_.PatternCount(), 2.0),
      using_likelihoods_(using_likelihoods),
      parsimony_pvs_(mmap_parsimony_path, 0, site_pattern_.PatternCount(), 2.0),
      parsimony_cost_matrix_(SankoffMatrix()),
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

  // Set Branch Length Optimization Method.
  optimization_method_ =
      (use_gradients ? Optimization::OptimizationMethod::BrentOptimizationWithGradients
                     : Optimization::OptimizationMethod::BrentOptimization);
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

// ** Scoring by Likelihood

double TPEngine::GetTopTreeLikelihoodWithEdge(const EdgeId edge_id) {
  return top_tree_log_likelihoods_per_edge_[edge_id.value_];
}

void TPEngine::InitializeLikelihood() {
  auto &pvs = likelihood_pvs_;
  // Set all PVs to Zero
  for (EdgeId edge_id = 0; edge_id < pvs.GetCount(); edge_id++) {
    for (const auto pv_type : PLVTypeEnum::Iterator()) {
      pvs.GetPV(pv_type, edge_id).setZero();
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
  const auto leafward_node_ids = dag_.LeafwardNodeTraversalTrace(true);
  for (const auto node_id : leafward_node_ids) {
    PopulateLeafwardLikelihoodPVForNode(node_id);
  }
}

void TPEngine::UpdateLikelihoodsAfterDAGAddNodePair(const NNIOperation &post_nni,
                                                    const NNIOperation &pre_nni,
                                                    std::optional<size_t> new_tree_id) {
  // Copy over branch lengths.
  CopyOverEdgeDataFromPreNNIToPostNNI(post_nni, pre_nni, new_tree_id);
  // Populate PVs.
  InitializeLikelihood();
}

double TPEngine::GetTopTreeLikelihoodWithProposedNNI(const NNIOperation &post_nni,
                                                     const NNIOperation &pre_nni,
                                                     const size_t spare_offset) {
  using NNIClade = NNIOperation::NNIClade;
  using NNICladeEnum = NNIOperation::NNICladeEnum;
  auto &pvs = likelihood_pvs_;
  GrowEdgeData(dag_.EdgeCountWithLeafSubsplits());
  // Node ids from pre-NNI in DAG.
  NNICladeEnum::Array<PVId> pre_id_map;
  NNICladeEnum::Array<PVId> post_id_map;
  const auto pre_parent_id = dag_.GetDAGNodeId(pre_nni.GetParent());
  const auto pre_child_id = dag_.GetDAGNodeId(pre_nni.GetChild());
  const auto pre_edge_id = dag_.GetEdgeIdx(pre_parent_id, pre_child_id);
  const auto pre_sister_clade =
      Bitset::Opposite(dag_.GetDAGEdge(pre_edge_id).GetSubsplitClade());
  // Create mapping between pre-NNI and post-NNI.
  const auto clade_map =
      NNIOperation::BuildNNICladeMapFromPreNNIToNNI(pre_nni, post_nni);
  // PLV ids from pre-NNI in DAG.
  auto choices = choice_map_.GetEdgeChoice(pre_edge_id);
  pre_id_map[NNIClade::ParentFocal] =
      pvs.GetPVIndex(PLVType::RHat, choices.parent_edge_id);
  pre_id_map[NNIClade::ParentSister] =
      pvs.GetPVIndex(PLVTypeEnum::PPLVType(pre_sister_clade), choices.parent_edge_id);
  pre_id_map[NNIClade::ChildLeft] = pvs.GetPVIndex(PLVType::PHatLeft, pre_edge_id);
  pre_id_map[NNIClade::ChildRight] = pvs.GetPVIndex(PLVType::PHatRight, pre_edge_id);
  // Use clade mapping to reference pre-NNI PVs for post-NNI PVs.
  for (const auto nni_clade : NNICladeEnum::Iterator()) {
    post_id_map[nni_clade] = pre_id_map[clade_map[nni_clade]];
  }
  // Get temp locations for post-NNI PVs and edge lengths.
  auto child_phat_pvid = pvs.GetSparePVIndex((spare_offset * 2));
  auto parent_rfocal_pvid = pvs.GetSparePVIndex((spare_offset * 2) + 1);
  auto post_edge_id = EdgeId(GetEdgeCount() + spare_offset);
  branch_lengths_[post_edge_id.value_] = branch_lengths_[pre_edge_id.value_];
  // Evolve post-parent with post-sister together.
  MultiplyPVs(parent_rfocal_pvid, post_id_map[NNIClade::ParentFocal],
              post_id_map[NNIClade::ParentSister]);
  // Evolve post-leftchild with post-rightchild together.
  MultiplyPVs(child_phat_pvid, post_id_map[NNIClade::ChildLeft],
              post_id_map[NNIClade::ChildRight]);
  // Get likelihood by evolving up from child to parent.
  ComputeLikelihood(post_edge_id, child_phat_pvid, parent_rfocal_pvid);
  double top_tree_log_likelihood = GetTopTreeLikelihoodWithEdge(post_edge_id);
  return top_tree_log_likelihood;
}

// For rootward traversal. Compute the likelihood PV for a given node, by using
// the left and right child edges from the choice map of the edge below node.
void TPEngine::PopulateRootwardLikelihoodPVForNode(const NodeId node_id) {
  // Populate edge PLV by evolving up given edge.
  auto PopulateEdgePLV = [this](const EdgeId edge_id) {
    auto &pvs = likelihood_pvs_;
    const auto choices = choice_map_.GetEdgeChoice(edge_id);
    // Evolve up from left child.
    if (choices.left_child_edge_id != NoId) {
      EvolveLikelihoodPPVUpEdge(edge_id, choices.left_child_edge_id);
    }
    // Evolve up from right child.
    if (choices.right_child_edge_id != NoId) {
      EvolveLikelihoodPPVUpEdge(edge_id, choices.right_child_edge_id);
    }
    // Update P-PLV from PHatLeft and PHatRight.
    const PVId p_pvid = pvs.GetPVIndex(PLVType::P, edge_id);
    const PVId phatleft_pvid = pvs.GetPVIndex(PLVType::PHatLeft, edge_id);
    const PVId phatright_pvid = pvs.GetPVIndex(PLVType::PHatRight, edge_id);
    if ((choices.left_child_edge_id != NoId) && (choices.right_child_edge_id != NoId)) {
      MultiplyPVs(p_pvid, phatleft_pvid, phatright_pvid);
    } else if (choices.right_child_edge_id != NoId) {
      TakePVValue(p_pvid, phatright_pvid);
    } else if (choices.left_child_edge_id != NoId) {
      TakePVValue(p_pvid, phatleft_pvid);
    }
  };

  for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto adj_node_id :
         dag_.GetDAGNode(node_id).GetNeighbors(Direction::Rootward, clade)) {
      const EdgeId edge_id = dag_.GetEdgeIdx(adj_node_id, node_id);
      PopulateEdgePLV(edge_id);
    }
  }
}

// For leafward traversal. Compute the likelihood PV for a given node, by using the
// parent and sister edges from the choice map of the edge above node.
void TPEngine::PopulateLeafwardLikelihoodPVForNode(const NodeId node_id) {
  // Populate edge PLV by evolving down given edge.
  auto PopulateEdgePLV = [this](const EdgeId edge_id) {
    auto &pvs = likelihood_pvs_;
    const auto choices = choice_map_.GetEdgeChoice(edge_id);
    // Evolve down parent.
    if (choices.parent_edge_id != NoId) {
      EvolveLikelihoodRPVDownEdge(choices.parent_edge_id, edge_id);
    }
    // Evolve with sister.
    const PVId center_pvid = pvs.GetPVIndex(PLVType::RHat, edge_id);
    const PVId left_pvid = pvs.GetPVIndex(PLVType::RLeft, edge_id);
    const PVId right_pvid = pvs.GetPVIndex(PLVType::RRight, edge_id);
    const PVId left_child_pvid = pvs.GetPVIndex(PLVType::PHatLeft, edge_id);
    const PVId right_child_pvid = pvs.GetPVIndex(PLVType::PHatRight, edge_id);
    MultiplyPVs(left_pvid, right_child_pvid, center_pvid);
    MultiplyPVs(right_pvid, left_child_pvid, center_pvid);
  };

  for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto adj_node_id :
         dag_.GetDAGNode(node_id).GetNeighbors(Direction::Leafward, clade)) {
      const EdgeId edge_id = dag_.GetEdgeIdx(node_id, adj_node_id);
      PopulateEdgePLV(edge_id);
    }
  }
}

void TPEngine::PopulateLeafLikelihoodPVsWithSitePatterns() {
  auto &pvs = likelihood_pvs_;
  for (PVId pv_id = 0; pv_id < pvs.GetPVCount(); pv_id++) {
    pvs.GetPV(pv_id).setZero();
  }
  NodeId node_id = 0;
  for (const auto &pattern : site_pattern_.GetPatterns()) {
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative symbol!");
      for (const auto pv_type : {PLVType::P}) {
        const auto node = dag_.GetDAGNode(node_id);
        for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
          for (const auto adj_node_id : node.GetNeighbors(Direction::Rootward, clade)) {
            const auto edge_id = dag_.GetEdgeIdx(node_id, adj_node_id);
            if (symbol == MmappedNucleotidePLV::base_count_) {  // Gap character.
              pvs.GetPV(pvs.GetPVIndex(pv_type, edge_id)).col(site_idx).setConstant(1.);
            } else if (symbol < MmappedNucleotidePLV::base_count_) {
              pvs.GetPV(pvs.GetPVIndex(pv_type, edge_id))(symbol, site_idx) = 1.;
            }
          }
        }
      }
      site_idx++;
    }
    node_id++;
  }
}

void TPEngine::PopulateRootLikelihoodPVsWithStationaryDistribution() {
  const NodeId root_node_id = dag_.GetDAGRootNodeId();
  for (const auto pv_type : {PLVType::RHat, PLVType::RLeft, PLVType::RRight}) {
    auto &pvs = likelihood_pvs_;
    for (const auto &node_id : dag_.GetRootsplitNodeIds()) {
      const auto edge_id = dag_.GetEdgeIdx(root_node_id, node_id);
      auto &pv = pvs.GetPV(pv_type, edge_id);
      for (Eigen::Index row_idx = 0; row_idx < pv.rows(); ++row_idx) {
        pv.row(row_idx).array() = stationary_distribution_(row_idx);
      }
    }
  }
}

void TPEngine::ComputeLikelihoods() {
  auto &pvs = likelihood_pvs_;

  for (EdgeId edge_id = 0; edge_id < dag_.EdgeCountWithLeafSubsplits(); edge_id++) {
    const auto choices = choice_map_.GetEdgeChoice(edge_id);
    const auto edge = dag_.GetDAGEdge(edge_id);
    if (choices.parent_edge_id != NoId) {
      const PVId parent_pvid = pvs.GetPVIndex(
          PLVTypeEnum::RPLVType(edge.GetSubsplitClade()), choices.parent_edge_id);
      const PVId child_pvid = pvs.GetPVIndex(PLVType::P, edge_id);
      ComputeLikelihood(edge_id, child_pvid, parent_pvid);
    }
  }

  top_tree_log_likelihoods_per_edge_ =
      log_likelihoods_.block(0, 0, dag_.EdgeCountWithLeafSubsplits(),
                             log_likelihoods_.cols()) *
      site_pattern_weights_;
}

void TPEngine::EvolveLikelihoodPPVUpEdge(const EdgeId rootward_edge_id,
                                         const EdgeId current_edge_id) {
  auto &pvs = likelihood_pvs_;
  const auto edge = dag_.GetDAGEdge(current_edge_id);
  const PVId parent_pvid =
      pvs.GetPVIndex(PLVTypeEnum::PPLVType(edge.GetSubsplitClade()), rootward_edge_id);
  const PVId child_pvid = pvs.GetPVIndex(PLVType::P, current_edge_id);
  SetToEvolvedPV(parent_pvid, current_edge_id, child_pvid);
}

void TPEngine::EvolveLikelihoodRPVDownEdge(const EdgeId current_edge_id,
                                           const EdgeId leafward_edge_id) {
  auto &pvs = likelihood_pvs_;

  const auto edge = dag_.GetDAGEdge(leafward_edge_id);
  const PLVType parent_focal_rpv = PLVTypeEnum::RPLVType(edge.GetSubsplitClade());
  const PVId parent_pvid = pvs.GetPVIndex(parent_focal_rpv, current_edge_id);
  const PVId child_pvid = pvs.GetPVIndex(PLVType::RHat, leafward_edge_id);
  SetToEvolvedPV(child_pvid, leafward_edge_id, parent_pvid);
}

// ** Scoring by Parsimony

double TPEngine::GetTopTreeParsimonyWithEdge(const EdgeId edge_id) {
  return top_tree_parsimony_per_edge_[edge_id.value_];
}

void TPEngine::PopulateLeafParsimonyPVsWithSitePatterns() {
  // first check that the psv_handler has been resized to deal with the leaf labels
  Assert(parsimony_pvs_.GetCount() >= site_pattern_.TaxonCount(),
         "Error in SankoffHandler::GenerateLeafPartials: "
         "parsimony_pvs_ should be initialized to accomodate"
         "the number of leaf nodes in the site_pattern_.");

  // Iterate over all leaf nodes to instantiate each with P partial values
  for (const auto node_id : dag_.GetLeafNodeIds()) {
    for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
      for (const auto adj_node_id :
           dag_.GetDAGNode(node_id).GetNeighbors(Direction::Rootward, clade)) {
        const auto edge_id = dag_.GetEdgeIdx(adj_node_id, node_id);
        SankoffPartial leaf_partials(state_count_, site_pattern_.PatternCount());
        // set leaf node partial to have big_double_ infinity substitute
        leaf_partials.block(0, 0, state_count_, site_pattern_.PatternCount())
            .fill(big_double_);
        // now fill in appropriate entries of the leaf-partial where non-infinite
        for (size_t pattern_idx = 0; pattern_idx < site_pattern_.PatternCount();
             pattern_idx++) {
          auto site_val = site_pattern_.GetPatternSymbol(node_id.value_, pattern_idx);
          if (site_val < state_count_) {
            leaf_partials(site_val, pattern_idx) = 0.;
          } else if (site_val == state_count_) {
            // Leaves with gaps in sequence and ambiguous nucleotides are assigned
            // sankoff partial vector [0, 0, 0, 0] at the corresponding site.
            leaf_partials.col(pattern_idx).fill(0);
          } else {
            Failwith(
                "Error in SankoffHandler::GenerateLeafPartials: Invalid nucleotide "
                "state in sequence alignment.");
          }
        }
        parsimony_pvs_.GetPV(PSVType::PLeft, edge_id) = leaf_partials;
        parsimony_pvs_.GetPV(PSVType::PRight, edge_id).fill(0);
      }
    }
  }
}

void TPEngine::InitializeParsimony() {
  auto &pvs = parsimony_pvs_;
  // Set all PVs to Zero
  for (EdgeId edge_id = 0; edge_id < pvs.GetCount(); edge_id++) {
    for (const auto pv_type : PSVTypeEnum::Iterator()) {
      pvs.GetPV(pv_type, edge_id).setZero();
    }
  }
  // Populate Leaves with Site Patterns.
  PopulateLeafParsimonyPVsWithSitePatterns();

  // Rootward Pass (populate P PVs)
  const auto rootward_node_ids = dag_.RootwardNodeTraversalTrace(false);
  for (const auto node_id : rootward_node_ids) {
    PopulateRootwardParsimonyPVForNode(node_id);
  }
  // Leafward Pass(populate R PVs)
  const auto leafward_node_ids = dag_.LeafwardNodeTraversalTrace(true);
  for (const auto node_id : leafward_node_ids) {
    PopulateLeafwardParsimonyPVForNode(node_id);
  }
}

void TPEngine::UpdateParsimoniesAfterDAGAddNodePair(const NNIOperation &post_nni,
                                                    const NNIOperation &pre_nni,
                                                    std::optional<size_t> new_tree_id) {
  // Copy over branch lengths.
  CopyOverEdgeDataFromPreNNIToPostNNI(post_nni, pre_nni, new_tree_id);
  // Populate PVs.
  InitializeParsimony();
}

void TPEngine::PopulateRootwardParsimonyPVForNode(const NodeId node_id) {
  // Populate edge PLV by accumulating parsimony up given edge.
  auto PopulateEdgePSV = [this](const EdgeId edge_id) {
    const auto choices = choice_map_.GetEdgeChoice(edge_id);
    if (choices.left_child_edge_id != NoId && choices.right_child_edge_id != NoId) {
      // Accumulate parsimony from left and right child.
      const EdgeId left_child_edge_id = choices.left_child_edge_id;
      const EdgeId right_child_edge_id = choices.right_child_edge_id;

      PopulateRootwardParsimonyPVForEdge(edge_id, left_child_edge_id,
                                         right_child_edge_id);
    }
  };

  for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto adj_node_id :
         dag_.GetDAGNode(node_id).GetNeighbors(Direction::Rootward, clade)) {
      const EdgeId edge_id = dag_.GetEdgeIdx(adj_node_id, node_id);
      PopulateEdgePSV(edge_id);
    }
  }
}

void TPEngine::PopulateLeafwardParsimonyPVForNode(const NodeId node_id) {
  // Populate edge PLV by accumulating parsimony down given edge.
  auto PopulateEdgePSV = [this](const EdgeId edge_id) {
    const auto choices = choice_map_.GetEdgeChoice(edge_id);
    if (choices.parent_edge_id != NoId && choices.sister_edge_id != NoId) {
      // Evolve down parent.
      const auto edge = dag_.GetDAGEdge(edge_id);
      const EdgeId parent_edge_id = choices.parent_edge_id;
      const EdgeId sister_edge_id = choices.sister_edge_id;
      const EdgeId left_child_edge_id =
          (edge.GetSubsplitClade() == SubsplitClade::Left) ? edge_id : sister_edge_id;
      const EdgeId right_child_edge_id =
          (edge.GetSubsplitClade() == SubsplitClade::Left) ? sister_edge_id : edge_id;

      PopulateLeafwardParsimonyPVForEdge(parent_edge_id, left_child_edge_id,
                                         right_child_edge_id);
    }
  };

  for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto adj_node_id :
         dag_.GetDAGNode(node_id).GetNeighbors(Direction::Leafward, clade)) {
      const EdgeId edge_id = dag_.GetEdgeIdx(node_id, adj_node_id);
      PopulateEdgePSV(edge_id);
    }
  }
}

void TPEngine::ComputeParsimonies() {
  for (EdgeId edge_id = 0; edge_id < dag_.EdgeCountWithLeafSubsplits(); edge_id++) {
    top_tree_parsimony_per_edge_[edge_id.value_] = ParsimonyScore(edge_id);
  }
}

EigenVectorXd TPEngine::ParentPartial(EigenVectorXd child_partials) {
  Assert(child_partials.size() == state_count_,
         "child_partials in SankoffHandler::ParentPartial should have 4 states.");
  EigenVectorXd parent_partials(state_count_);
  parent_partials.setZero();
  for (size_t parent_state = 0; parent_state < state_count_; parent_state++) {
    EigenVectorXd partials_for_state(state_count_);
    for (size_t child_state = 0; child_state < state_count_; child_state++) {
      auto cost = parsimony_cost_matrix_.GetCost(parent_state, child_state);
      partials_for_state[child_state] = cost + child_partials[child_state];
    }
    auto minimum_element =
        *std::min_element(partials_for_state.data(),
                          partials_for_state.data() + partials_for_state.size());
    parent_partials[parent_state] = minimum_element;
  }

  return parent_partials;
}

EigenVectorXd TPEngine::TotalPPartial(EdgeId edge_id, size_t site_idx) {
  return parsimony_pvs_.GetPV(PSVType::PLeft, edge_id).col(site_idx) +
         parsimony_pvs_.GetPV(PSVType::PRight, edge_id).col(site_idx);
}

void TPEngine::PopulateRootwardParsimonyPVForEdge(const EdgeId parent_id,
                                                  const EdgeId left_child_id,
                                                  const EdgeId right_child_id) {
  for (size_t pattern_idx = 0; pattern_idx < site_pattern_.PatternCount();
       pattern_idx++) {
    // Which child partial is in right or left doesn't actually matter because they
    // are summed when calculating q_partials.
    parsimony_pvs_.GetPV(PSVType::PLeft, parent_id).col(pattern_idx) =
        ParentPartial(TotalPPartial(left_child_id, pattern_idx));

    parsimony_pvs_.GetPV(PSVType::PRight, parent_id).col(pattern_idx) =
        ParentPartial(TotalPPartial(right_child_id, pattern_idx));
  }
}

void TPEngine::PopulateLeafwardParsimonyPVForEdge(const EdgeId parent_id,
                                                  const EdgeId left_child_id,
                                                  const EdgeId right_child_id) {
  for (size_t pattern_idx = 0; pattern_idx < site_pattern_.PatternCount();
       pattern_idx++) {
    auto partials_from_parent =
        ParentPartial(parsimony_pvs_.GetPV(PSVType::Q, parent_id).col(pattern_idx));
    for (const auto child_id : {left_child_id, right_child_id}) {
      EdgeId sister_id = ((child_id == left_child_id) ? right_child_id : left_child_id);
      auto partials_from_sister = ParentPartial(TotalPPartial(sister_id, pattern_idx));
      parsimony_pvs_.GetPV(PSVType::Q, child_id).col(pattern_idx) =
          partials_from_sister + partials_from_parent;
    }
  }
}

double TPEngine::ParsimonyScore(EdgeId edge_id) {
  auto weights = site_pattern_.GetWeights();
  double total_parsimony = 0.;
  for (size_t pattern = 0; pattern < site_pattern_.PatternCount(); pattern++) {
    // Note: doing ParentPartial first for the left and right p_partials and then adding
    // them together will give the same minimum parsimony score, but doesn't give
    // correct Sankoff Partial vector for the new rooting
    auto total_tree = ParentPartial(TotalPPartial(edge_id, pattern));
    total_tree += ParentPartial(parsimony_pvs_.GetPV(PSVType::Q, edge_id).col(pattern));

    // If node_id is the root node, calculating the total_tree vector like so does not
    // yield the SankoffPartial of an actual rooting, but this will not change the
    // minimum value in the partial, so the root node can still be used to calculate the
    // parsimony score.
    total_parsimony +=
        *std::min_element(total_tree.begin(), total_tree.end()) * weights[pattern];
  }
  return total_parsimony;
}

// ** PV Operations for Likelihoods

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

// ** Branch Length Optimization

void TPEngine::SetOptimizationMethod(const Optimization::OptimizationMethod method) {
  optimization_method_ = method;
}

void TPEngine::Optimization(const EdgeId edge_id) {
  switch (optimization_method_) {
    case Optimization::OptimizationMethod::BrentOptimization:
      return BrentOptimization(op);
    case Optimization::OptimizationMethod::BrentOptimizationWithGradients:
      return BrentOptimizationWithGradients(op);
    case Optimization::OptimizationMethod::GradientAscentOptimization:
      return GradientAscentOptimization(op);
    case Optimization::OptimizationMethod::LogSpaceGradientAscentOptimization:
      return LogSpaceGradientAscentOptimization(op);
    case Optimization::OptimizationMethod::NewtonOptimization:
      return NewtonOptimization(op);
    default:
      Failwith("TPEngine::Optimization(): Invalid OptimizationMethod given.");
  }
}

void TPEngine::SetSignificantDigitsForOptimization(int significant_digits) {
  significant_digits_for_optimization_ = significant_digits;
}

void TPEngine::BrentOptimization(const EdgeId edge_id) {
  if (branch_length_differences_(edge_id.value_) <
      branch_length_difference_threshold_) {
    return;
  }

  auto edge = dag_.GetDAGEdge(edge_id);

  auto negative_log_likelihood = [this, edge](double log_branch_length) {
    SetTransitionMatrixToHaveBranchLength(exp(log_branch_length));
    PreparePerPatternLogLikelihoodsForEdge(edge.GetParent(), edge.GetChild());
    return -per_pattern_log_likelihoods_.dot(site_pattern_weights_);
  };

  double current_log_branch_length = log(branch_lengths_(edge_id.value_));
  double current_neg_log_likelihood =
      negative_log_likelihood(current_log_branch_length);

  const auto [log_branch_length, neg_log_likelihood] =
      Optimization::BrentMinimize<false>(
          negative_log_likelihood, current_log_branch_length, min_log_branch_length_,
          max_log_branch_length_, significant_digits_for_optimization_,
          max_iter_for_optimization_, step_size_for_log_space_optimization_);

  // Numerical optimization sometimes yields new nllk > current nllk.
  // In this case, we reset the branch length to the previous value.
  if (neg_log_likelihood > current_neg_log_likelihood) {
    branch_lengths_(edge_id.value_) = exp(current_log_branch_length);
  } else {
    branch_lengths_(edge_id.value_) = exp(log_branch_length);
  }
  branch_length_differences_(edge_id.value_) =
      abs(exp(current_log_branch_length) - branch_lengths_(edge_id.value_));
}

void TPEngine::BrentOptimizationWithGradients(const EdgeId edge_id) {
  if (branch_length_differences_(edge_id.value_) <
      branch_length_difference_threshold_) {
    return;
  }

  auto negative_log_likelihood_and_derivative = [this,
                                                 edge_id](double log_branch_length) {
    double branch_length = exp(log_branch_length);
    branch_lengths_(edge_id.value_) = branch_length;
    auto [log_likelihood, log_likelihood_derivative] =
        this->LogLikelihoodAndDerivative(edge_id);
    return std::make_pair(-log_likelihood, -branch_length * log_likelihood_derivative);
  };

  double current_log_branch_length = log(branch_lengths_(edge_id.value_));
  double current_neg_log_likelihood =
      negative_log_likelihood_and_derivative(current_log_branch_length).first;
  const auto [log_branch_length, neg_log_likelihood] =
      Optimization::BrentMinimize<true>(
          negative_log_likelihood_and_derivative, current_log_branch_length,
          min_log_branch_length_, max_log_branch_length_,
          significant_digits_for_optimization_, max_iter_for_optimization_,
          step_size_for_log_space_optimization_);

  if (neg_log_likelihood > current_neg_log_likelihood) {
    branch_lengths_(edge_id.value_) = exp(current_log_branch_length);
  } else {
    branch_lengths_(edge_id.value_) = exp(log_branch_length);
  }
  branch_length_differences_(edge_id.value_) =
      abs(exp(current_log_branch_length) - branch_lengths_(edge_id.value_));
}

void TPEngine::GradientAscentOptimization(const EdgeId edge_id) {
  auto log_likelihood_and_derivative = [this, &op](double branch_length) {
    branch_lengths_(op.gpcsp_) = branch_length;
    return this->LogLikelihoodAndDerivative(op);
  };
  const auto branch_length = Optimization::GradientAscent(
      log_likelihood_and_derivative, branch_lengths_(op.gpcsp_),
      significant_digits_for_optimization_, step_size_for_optimization_,
      min_log_branch_length_, max_iter_for_optimization_);
  branch_lengths_(op.gpcsp_) = branch_length;
}

void TPEngine::LogSpaceGradientAscentOptimization(const EdgeId edge_id) {
  auto log_likelihood_and_derivative = [this, edge_id](double branch_length) {
    branch_lengths_(edge_id.value_) = branch_length;
    return this->LogLikelihoodAndDerivative(edge_id);
  };
  const auto branch_length = Optimization::LogSpaceGradientAscent(
      log_likelihood_and_derivative, branch_lengths_(edge_id.value_),
      significant_digits_for_optimization_, step_size_for_log_space_optimization_,
      exp(min_log_branch_length_), max_iter_for_optimization_);
  branch_lengths_(edge_id.value_) = branch_length;
}

void TPEngine::NewtonOptimization(const EdgeId edge_id) {
  if (branch_length_differences_(edge_id.value_) <
      branch_length_difference_threshold_) {
    return;
  }

  auto log_likelihood_and_first_two_derivatives = [this,
                                                   edge_id](double log_branch_length) {
    double x = exp(log_branch_length);
    branch_lengths_(edge_id.value_) = x;
    auto [f_x, f_prime_x, f_double_prime_x] =
        this->LogLikelihoodAndFirstTwoDerivatives(edge_id);
    // x = exp(y) --> f'(exp(y)) = exp(y) * f'(exp(y)) = x * f'(x)
    double f_prime_y = x * f_prime_x;
    double f_double_prime_y = f_prime_y + std::pow(x, 2) * f_double_prime_x;
    return std::make_tuple(f_x, f_prime_y, f_double_prime_y);
  };

  double current_log_branch_length = log(branch_lengths_(edge_id.value_));
  const auto log_branch_length = Optimization::NewtonRaphsonOptimization(
      log_likelihood_and_first_two_derivatives, current_log_branch_length,
      significant_digits_for_optimization_, denominator_tolerance_for_newton_,
      min_log_branch_length_, max_log_branch_length_, max_iter_for_optimization_);

  branch_lengths_(edge_id.value_) = exp(log_branch_length);

  branch_length_differences_(edge_id.value_) =
      abs(exp(current_log_branch_length) - branch_lengths_(edge_id.value_));
}

// ** Parameter Data

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
  }
  // Reindex work space to realign with DAG.
  if (node_reindexer.has_value()) {
    ReindexNodeData(node_reindexer.value(), old_node_count);
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
    if (using_likelihoods_) {
      likelihood_pvs_.Resize(new_edge_count, GetAllocatedEdgeCount());
    }
    if (using_parsimony_) {
      parsimony_pvs_.Resize(new_edge_count, GetAllocatedEdgeCount());
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
         "Node Reindexer is the wrong size for TPEngine.");
  Assert(node_reindexer.IsValid(GetNodeCount()), "Node Reindexer is not valid.");
}

void TPEngine::ReindexEdgeData(const Reindexer &edge_reindexer,
                               const size_t old_edge_count) {
  Assert(edge_reindexer.size() == GetEdgeCount(),
         "Edge Reindexer is the wrong size for TPEngine.");
  Assert(edge_reindexer.IsValid(GetEdgeCount()),
         "Edge Reindexer is not valid for TPEngine size.");
  // Reindex data vectors.
  Reindexer::ReindexInPlace<EigenVectorXd, double>(branch_lengths_, edge_reindexer,
                                                   GetEdgeCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(top_tree_log_likelihoods_per_edge_,
                                                   edge_reindexer, GetEdgeCount());
  // Reindex PVs.
  Reindexer pv_likelihood_reindexer =
      likelihood_pvs_.BuildPVReindexer(edge_reindexer, old_edge_count, GetEdgeCount());
  likelihood_pvs_.Reindex(pv_likelihood_reindexer);
  Reindexer pv_parsimony_reindexer =
      likelihood_pvs_.BuildPVReindexer(edge_reindexer, old_edge_count, GetEdgeCount());
  parsimony_pvs_.Reindex(pv_parsimony_reindexer);
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
