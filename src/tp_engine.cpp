// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "tp_engine.hpp"
#include "gp_engine.hpp"

TPEngine::TPEngine(GPDAG &dag, SitePattern &site_pattern,
                   const std::string &mmap_file_path, bool using_likelihoods,
                   bool using_parsimony)
    : dag_(dag),
      site_pattern_(site_pattern),
      likelihood_pvs_(mmap_file_path, 0, site_pattern_.PatternCount(), 2.0),
      using_likelihoods_(using_likelihoods),
      parsimony_pvs_(mmap_file_path, 0, site_pattern_.PatternCount(), 2.0),
      using_parsimony_(using_parsimony),
      choice_map_(dag) {
  // Initialize site pattern-based data.
  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ = EigenVectorXdOfStdVectorDouble(weights);
  // Initialize node-based data
  GrowNodeData(dag_.NodeCountWithoutDAGRoot(), std::nullopt, std::nullopt, true);
  // Initialize edge-based data
  GrowEdgeData(dag_.EdgeCountWithLeafSubsplits(), std::nullopt, std::nullopt, true);
  // Initialize scores.
  InitializeChoiceMap();
}

// ** Initialization

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

void TPEngine::ComputeTopTreeLikelihoodWithEdge(const EdgeId edge_id) {
  // const auto edge = dag_.GetDAGEdge(edge_id);
  // SetTransitionMatrixToHaveBranchLength(branch_lengths_(edge_id.value_));
  // PreparePerPatternLogLikelihoodsForEdge(edge.GetParent().value_,
  //                                        edge.GetChild().value_);
  // log_likelihoods_.row(edge_id.value_) = per_pattern_log_likelihoods_;
}

double TPEngine::GetTopTreeLikelihoodWithEdge(const EdgeId edge_id) {
  return top_tree_log_likelihoods_per_edge_[edge_id.value_];
}

void TPEngine::InitializeLikelihood() {
  std::cout << "INIT_LIKELIHOOD [BEGIN]" << std::endl;
  auto &pvs = likelihood_pvs_;
  // Set all PVs to Zero
  for (NodeId node_id = 0; node_id < pvs.GetNodeCount(); node_id++) {
    for (const auto pv_type : PLVTypeEnum::Iterator()) {
      pvs.GetPV(PLVType::RHat, node_id).setZero();
    }
  }
  // Populate Leaves with Site Patterns.
  PopulateLeafLikelihoodPVsWithSitePatterns();
  // Populate Rootsplit with Stationary Distribution.
  PopulateRootLikelihoodPVsWithStationaryDistribution();
  // Rootward Pass (populate P PVs)
  const auto rootward_node_ids = dag_.RootwardNodeTraversalTrace(false);
  std::cout << "rootward_trace: " << rootward_node_ids << std::endl;
  for (const auto node_id : rootward_node_ids) {
    PopulateRootwardLikelihoodPVForNode(node_id);
  }
  // Leafward Pass (populate R PVs)
  const auto leafward_node_ids = dag_.LeafwardNodeTraversalTrace(false);
  for (const auto node_id : leafward_node_ids) {
    PopulateLeafwardLikelihoodPVForNode(node_id);
  }

  std::cout << "INIT_LIKELIHOOD [END]" << std::endl;
}

void TPEngine::UpdateDAGAfterAddNodePairByLikelihood(const NNIOperation &nni_op) {
  // !ISSUE #440
  Failwith("Currently no implementation.");
}

// For rootward traversal. Compute the likelihood PV for a given node, by using
// the left and right child edges from the choice map of the edge below node.
void TPEngine::PopulateRootwardLikelihoodPVForNode(const NodeId node_id) {
  EdgeIdVector central_edge_ids;
  const auto node = dag_.GetDAGNode(node_id);
  // Iterate over all edges where the given node is the child_node of the edge.
  size_t option_count = 0;
  for (const auto focal_clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto rootward_node_id :
         node.GetNeighbors(Direction::Rootward, focal_clade)) {
      option_count++;

      const auto edge_id = dag_.GetEdgeIdx(node_id, rootward_node_id);
      const auto edge = dag_.GetDAGEdge(edge_id);
      const auto edge_choice = choice_map_.GetEdgeChoice(edge_id);
      const PVId pv_center_index = likelihood_pvs_.GetPVIndex(PLVType::P, node_id);
      const PVId pv_left_index = likelihood_pvs_.GetPVIndex(PLVType::PHatLeft, node_id);
      const PVId pv_right_index =
          likelihood_pvs_.GetPVIndex(PLVType::PHatRight, node_id);

      std::cout << "ROOTWARD: " << node_id << std::endl << edge_choice << std::endl;

      // Evolve up from left child.
      const EdgeId left_child_edge_id = edge_choice.left_child_edge_id;
      if (left_child_edge_id != NoId) {
        const auto left_child_edge = dag_.GetDAGEdge(left_child_edge_id);
        const NodeId left_child_node_id = left_child_edge.GetChild();
        // left_child_p(s)
        const PVId pv_left_child_index =
            likelihood_pvs_.GetPVIndex(PLVType::P, left_child_node_id);
        SetToEvolvedPV(pv_left_index, left_child_edge_id, pv_left_child_index);
      }

      // Evolve up from right child.
      const auto right_child_edge_id = edge_choice.right_child_edge_id;
      if (right_child_edge_id != NoId) {
        const auto right_child_edge = dag_.GetDAGEdge(right_child_edge_id);
        const NodeId right_child_node_id = right_child_edge.GetChild();
        // right_child_p(s)
        const PVId pv_right_child_index =
            likelihood_pvs_.GetPVIndex(PLVType::P, right_child_node_id);
        SetToEvolvedPV(pv_right_index, right_child_edge_id, pv_right_child_index);
      }

      if (right_child_edge_id != NoId && left_child_edge_id != NoId) {
        Multiply(pv_center_index, pv_left_index, pv_right_index);
      } else if (right_child_edge_id != NoId) {
        Set(pv_center_index, pv_right_index);
      } else if (left_child_edge_id != NoId) {
        Set(pv_center_index, pv_left_index);
      } else {
        std::cout << "NO_MAPPING: " << node_id << std::endl;
      }
    }
  }
  std::cout << "option_count [leafward]: " << node_id << " " << option_count
            << std::endl;
}

// For leafward traversal. Compute the likelihood PV for a given node, by using the
// parent and sister edges from the choice map of the edge above node.
void TPEngine::PopulateLeafwardLikelihoodPVForNode(const NodeId node_id) {
  auto &pvs = likelihood_pvs_;
  EdgeIdVector central_edge_ids;
  const auto node = dag_.GetDAGNode(node_id);
  // Iterate over all edges where the given node is the parent_node of the edge.
  size_t option_count = 0;
  for (const auto focal_clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto rootward_node_id :
         node.GetNeighbors(Direction::Leafward, focal_clade)) {
      option_count++;

      const PLVType focal_pv =
          (focal_clade == SubsplitClade::Left) ? PLVType::RLeft : PLVType::RRight;
      const PLVType sister_pv =
          (focal_clade == SubsplitClade::Left) ? PLVType::RRight : PLVType::RLeft;
      const PVId pv_center_index =
          likelihood_pvs_.GetPVIndex(PLVType::RHat, rootward_node_id);
      const PVId pv_focal_index =
          likelihood_pvs_.GetPVIndex(focal_pv, rootward_node_id);
      const PVId pv_sister_index =
          likelihood_pvs_.GetPVIndex(sister_pv, rootward_node_id);
      const auto edge_id = dag_.GetEdgeIdx(node_id, rootward_node_id);
      const auto edge = dag_.GetDAGEdge(edge_id);
      const auto edge_choice = choice_map_.GetEdgeChoice(edge_id);

      // Evolve down parent.
      const auto parent_edge_id = edge_choice.parent_edge_id;
      if (parent_edge_id != NoId) {
        const auto parent_edge = dag_.GetDAGEdge(parent_edge_id);
        // parent_q(s)
        const PVId pv_parent_parent_index =
            likelihood_pvs_.GetPVIndex(PLVType::RHat, parent_edge.GetParent());
        const PVId pv_parent_child_index =
            likelihood_pvs_.GetPVIndex(PLVType::RHat, parent_edge.GetChild());
        SetToEvolvedPV(pv_center_index, parent_edge_id, pv_parent_parent_index);
      }

      // Evolve up from sister.
      const auto sister_edge_id = edge_choice.sister_edge_id;
      if (sister_edge_id != NoId) {
        const auto sister_edge = dag_.GetDAGEdge(sister_edge_id);
        const NodeId sister_node_id = sister_edge.GetChild();
        // sister_p(s)
        const PVId pv_sister_child_index =
            likelihood_pvs_.GetPVIndex(sister_pv, sister_node_id);
        SetToEvolvedPV(pv_sister_index, sister_edge_id, pv_sister_child_index);
        Multiply(pv_sister_index, pv_sister_index, pv_center_index);
      }
    }
  }
  // std::cout << "option_count [leafward]: " << node_id << " " << option_count
  //           << std::endl;
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
  for (const auto &node_id : dag_.GetRootsplitNodeIds()) {
    auto &pv = likelihood_pvs_.GetPV(PLVType::RHat, node_id);
    for (Eigen::Index row_idx = 0; row_idx < pv.rows(); ++row_idx) {
      pv.row(row_idx).array() = stationary_distribution_(row_idx);
    }
  }
}

void TPEngine::ComputeLikelihoods() {
  std::cout << "ComputeLikelihoods [BEGIN]" << std::endl;
  dag_.IterateOverRealNodes([this](SubsplitDAGNode node) {
    dag_.IterateOverLeafwardEdges(
        node, [this, node](const bool is_edge_on_left, SubsplitDAGNode child_node) {
          auto &pvs = likelihood_pvs_;
          const EdgeId edge_id = dag_.GetEdgeIdx(node.Id(), child_node.Id());
          const auto edge = dag_.GetDAGEdge(edge_id);
          const PVId parent_pv_index =
              pvs.GetPVIndex(PLVHandler::RPLVType(is_edge_on_left), node.Id());
          const PVId child_pv_index = pvs.GetPVIndex(PLVType::P, child_node.Id());
          ComputeLikelihood(edge_id, parent_pv_index, child_pv_index);
        });
  });

  top_tree_log_likelihoods_per_edge_ =
      log_likelihoods_.block(0, 0, dag_.EdgeCountWithLeafSubsplits(),
                             log_likelihoods_.cols()) *
      site_pattern_weights_;
  std::cout << "TP_LIKELIHOODS: " << top_tree_log_likelihoods_per_edge_ << std::endl;
  std::cout << "TP_SITE_PATTERN: " << site_pattern_weights_ << std::endl;
  std::cout << "ComputeLikelihoods [END]" << std::endl;
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

// ** Operations

void TPEngine::SetTransitionMatrixToHaveBranchLength(const double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void TPEngine::Set(const PVId dest_id, const PVId src_id) {
  auto &pvs = likelihood_pvs_;
  pvs.GetPV(dest_id).array() = pvs.GetPV(src_id).array();
}

void TPEngine::Multiply(const PVId dest_id, const PVId src1_id, const PVId src2_id) {
  auto &pvs = likelihood_pvs_;
  pvs.GetPV(dest_id).array() = pvs.GetPV(src1_id).array() * pvs.GetPV(src2_id).array();
  // rescaling_counts_(op.dest_) =
  //     rescaling_counts_(op.src1_) + rescaling_counts_(op.src2_);
  // AssertPLVIsFinite(op.dest_, "Multiply dest_ is not finite");
  // RescalePLVIfNeeded(op.dest_);
}

void TPEngine::ComputeLikelihood(const EdgeId dest_id, const PVId child_id,
                                 const PVId parent_id) {
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(dest_id.value_));
  PreparePerPatternLogLikelihoodsForEdge(parent_id, child_id);
  log_likelihoods_.row(dest_id.value_) = per_pattern_log_likelihoods_;
  // std::cout << "TP_LIKE: " << branch_lengths_(dest_id.value_) << " " << dest_id << "
  // "
  //           << parent_id << " " << child_id << " " << per_pattern_log_likelihoods_
  //           << " " << transition_matrix_.norm() << std::endl;
  // std::cout << "TP_PARENT_PV: " << parent_id << std::endl
  //           << likelihood_pvs_.ToString(parent_id) << std::endl;
  // std::cout << "TP_CHILD_PV: " << child_id << std::endl
  //           << likelihood_pvs_.ToString(child_id) << std::endl;
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

void TPEngine::MultiplyAndEvolvePV(const PVId dest_id, const EdgeId edge_id,
                                   const PVId src1_id, const PVId src2_id) {
  auto &pvs = likelihood_pvs_;
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(edge_id.value_));
  Multiply(dest_id, src1_id, src2_id);
  pvs.GetPV(dest_id).array() = (transition_matrix_ * pvs.GetPV(dest_id)).array();
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
  }
  // Resize to fit without deallocating unused memory.
  branch_lengths_.conservativeResize(GetPaddedEdgeCount());
  log_likelihoods_.conservativeResize(GetPaddedEdgeCount(),
                                      site_pattern_.PatternCount());
  top_tree_log_likelihoods_per_edge_.conservativeResize(GetPaddedEdgeCount());

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

void TPEngine::FunctionOverRootedTreeCollection(
    FunctionOnTreeNodeByEdge function_on_tree_node_by_edge,
    const RootedTreeCollection &tree_collection, const BitsetSizeMap &edge_indexer) {
  const auto leaf_count = tree_collection.TaxonCount();
  const size_t default_index = branch_lengths_.size();
  for (const auto &tree : tree_collection.Trees()) {
    tree.Topology()->RootedPCSPPreorder(
        [&leaf_count, &default_index, &edge_indexer, &tree,
         &function_on_tree_node_by_edge](
            const Node *sister_node, const Node *focal_node,
            const Node *left_child_node, const Node *right_child_node) {
          Bitset edge_bitset =
              SBNMaps::PCSPBitsetOf(leaf_count, sister_node, false, focal_node, false,
                                    left_child_node, false, right_child_node, false);
          const auto edge_idx = AtWithDefault(edge_indexer, edge_bitset, default_index);
          if (edge_idx != default_index) {
            function_on_tree_node_by_edge(EdgeId(edge_idx), tree, focal_node);
          }
        },
        true);
  }
}

void TPEngine::SetChoiceMapByTakingFirst(const RootedTreeCollection &tree_collection,
                                         const BitsetSizeMap &edge_indexer) {
  size_t unique_edge_count = GetEdgeCount();
  branch_lengths_.setZero();
  EigenVectorXi observed_edge_counts = EigenVectorXi::Zero(unique_edge_count);
  auto set_choice_map_and_increment_edge_count =
      [&observed_edge_counts, this](const EdgeId edge_idx, const RootedTree &tree,
                                    const Node *focal_node) {
        branch_lengths_(edge_idx.value_) += tree.BranchLength(focal_node);
        observed_edge_counts(edge_idx.value_)++;
      };
  FunctionOverRootedTreeCollection(set_choice_map_and_increment_edge_count,
                                   tree_collection, edge_indexer);
  for (EdgeId edge_idx = 0; edge_idx < unique_edge_count; edge_idx++) {
    if (observed_edge_counts(edge_idx.value_) == 0) {
      branch_lengths_(edge_idx.value_) = default_branch_length_;
    } else {
      // Normalize the branch length total using the counts to get a mean branch
      // length.
      branch_lengths_(edge_idx.value_) /=
          static_cast<double>(observed_edge_counts(edge_idx.value_));
    }
  }
}

void TPEngine::SetBranchLengthByTakingFirst(const RootedTreeCollection &tree_collection,
                                            const BitsetSizeMap &edge_indexer) {
  // Unique edges in collection should be the same as the number of total edges in DAG
  // created from collection.
  size_t unique_edge_count = GetEdgeCount();
  branch_lengths_.setZero();
  EigenVectorXi observed_edge_counts = EigenVectorXi::Zero(unique_edge_count);
  // Set branch lengths on first occurance.
  auto set_first_branch_length_and_increment_edge_count =
      [&observed_edge_counts, this](const EdgeId edge_idx, const RootedTree &tree,
                                    const Node *focal_node) {
        if (observed_edge_counts(edge_idx.value_) == 0) {
          branch_lengths_(edge_idx.value_) = tree.BranchLength(focal_node);
          observed_edge_counts(edge_idx.value_)++;
        }
      };
  FunctionOverRootedTreeCollection(set_first_branch_length_and_increment_edge_count,
                                   tree_collection, edge_indexer);
  // Either set branch length to the first occurance if edge exists in collection,
  // otherwise set length to default.
  for (EdgeId edge_idx = 0; edge_idx < unique_edge_count; edge_idx++) {
    if (observed_edge_counts(edge_idx.value_) == 0) {
      branch_lengths_(edge_idx.value_) = default_branch_length_;
    }
  }
}
