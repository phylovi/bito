// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "nni_engine.hpp"

using PLVType = PLVHandler::PLVType;

NNIEngine::NNIEngine(GPDAG &dag, GPEngine *gp_engine, TPEngine *tp_engine)
    : dag_(dag), graft_dag_(std::make_unique<GraftDAG>(dag)) {
  if (gp_engine) {
    gp_engine_ = gp_engine;
  }
  if (tp_engine) {
    tp_engine_ = tp_engine;
  }
}

// ** Runners

void NNIEngine::Run() {
  RunInit();
  // Loop until no more eligible NNIs are found.
  while (GetAdjacentNNICount() > 0) {
    RunMainLoop();
    RunPostLoop();
  }
}

void NNIEngine::RunInit() {
  // Initialize Adjacent NNIs based on starting state of DAG.
  ResetAllNNIs();
  SyncAdjacentNNIsWithDAG();
  PrepGPEngineForLikelihoods();
  FilterInit();
}

void NNIEngine::RunMainLoop() {
  // (1) Add all adjacent NNIs to the GraftDAG.
  GraftAdjacentNNIsToDAG();
  // (2) Evaluate each adjacent NNI.
  FilterPreUpdate();
  FilterEvaluateAdjacentNNIs();
  FilterPostUpdate();
  // (3) Select whether to accept or reject adjacent NNIs via filter.
  FilterProcessAdjacentNNIs();
  // (4) Add accepted NNIs to permanent DAG.
  AddAcceptedNNIsToDAG();

  sweep_count_++;
}

void NNIEngine::RunPostLoop() {
  // (5) Update Adjacent NNIs to reflect added NNI.
  UpdateAdjacentNNIs(false);
  // (6) Reset Accepted NNIs and GraftDAG and save results.
  UpdateAcceptedNNIs(true);
  UpdateRejectedNNIs(true);
}

// ** Filter Functions

void NNIEngine::SetNoEvaluate() {
  filter_eval_fn_ = [](NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                       TPEngine &this_tp_engine, GraftDAG &this_graft_dag,
                       const NNIOperation &nni) -> double { return 0.0; };
}

void NNIEngine::SetNoFilter(const bool set_nni_to_pass) {
  filter_process_fn_ = [set_nni_to_pass](
                           NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                           TPEngine &this_tp_engine, GraftDAG &this_graft_dag,
                           const NNIOperation &nni,
                           const double nni_score) -> bool { return set_nni_to_pass; };
}

void NNIEngine::SetScoreCutoff(const double score_cutoff) {
  filter_process_fn_ = [score_cutoff](
                           NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                           TPEngine &this_tp_engine, GraftDAG &this_graft_dag,
                           const NNIOperation &nni, const double nni_score) -> bool {
    bool nni_passes = (nni_score >= score_cutoff);
    return nni_passes;
  };
}

// ** Filtering

void NNIEngine::FilterInit() {
  if (filter_init_fn_) {
    (filter_init_fn_)(*this, GetGPEngine(), GetTPEngine(), GetGraftDAG());
  }
}

void NNIEngine::FilterPreUpdate() {
  if (filter_pre_update_fn_) {
    (filter_pre_update_fn_)(*this, GetGPEngine(), GetTPEngine(), GetGraftDAG());
  }
}

void NNIEngine::FilterEvaluateAdjacentNNIs() {
  if (!filter_eval_fn_) return;
  for (const auto &nni : GetAdjacentNNIs()) {
    const double nni_score =
        (filter_eval_fn_)(*this, GetGPEngine(), GetTPEngine(), GetGraftDAG(), nni);
    AddScoreForNNI(nni, nni_score);
  }
}

void NNIEngine::FilterPostUpdate() {
  if (filter_post_update_fn_) {
    (filter_post_update_fn_)(*this, GetGPEngine(), GetTPEngine(), GetGraftDAG());
  }
}

void NNIEngine::FilterProcessAdjacentNNIs() {
  Assert(filter_process_fn_, "Must assign a filter process function.");
  for (const auto &nni : GetAdjacentNNIs()) {
    double nni_score = (*GetScoredNNIs().find(nni)).second;
    const bool accept_nni = (filter_process_fn_)(*this, GetGPEngine(), GetTPEngine(),
                                                 GetGraftDAG(), nni, nni_score);
    if (accept_nni) {
      accepted_nnis_.insert(nni);
    } else {
      rejected_nnis_.insert(nni);
    }
  }
}

// ** Filtering Scheme

void NNIEngine::SetFilterInitFunction(StaticFilterInitFunction filter_init_fn) {
  filter_init_fn_ = filter_init_fn;
}
void NNIEngine::SetFilterPreUpdateFunction(
    StaticFilterUpdateFunction filter_pre_update_fn) {
  filter_pre_update_fn_ = filter_pre_update_fn;
}
void NNIEngine::SetFilterEvalFunction(StaticFilterEvaluateFunction filter_eval_fn) {
  filter_eval_fn_ = filter_eval_fn;
}
void NNIEngine::SetFilterPostUpdateFunction(
    StaticFilterUpdateFunction filter_post_update_fn) {
  filter_post_update_fn_ = filter_post_update_fn;
}
void NNIEngine::SetFilterProcessFunction(
    StaticFilterProcessFunction filter_process_fn) {
  filter_process_fn_ = filter_process_fn;
}

void NNIEngine::SetGPLikelihoodFilteringScheme(const double score_cutoff) {
  // Pre-Update Evaluation function
  filter_pre_update_fn_ = [](NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                             TPEngine &this_tp_engine, GraftDAG &this_graft_dag) {
    this_nni_engine.ScoreAdjacentNNIsByGPLikelihood();
  };
  SetScoreCutoff(score_cutoff);
}

void NNIEngine::SetTPLikelihoodFilteringScheme(const double score_cutoff) {
  // Pre-Update Evaluation function
  filter_pre_update_fn_ = [](NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                             TPEngine &this_tp_engine, GraftDAG &this_graft_dag) {
    this_nni_engine.ScoreAdjacentNNIsByTPLikelihood();
  };
  SetScoreCutoff(score_cutoff);
}

void NNIEngine::SetTPParsimonyFilteringScheme(const double score_cutoff) {
  // Pre-Update Evaluation function
  filter_pre_update_fn_ = [](NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                             TPEngine &this_tp_engine, GraftDAG &this_graft_dag) {
    this_nni_engine.ScoreAdjacentNNIsByTPParsimony();
  };
  SetScoreCutoff(score_cutoff);
}

// ** Key Indexing

NNIEngine::KeyIndex NNIEngine::NNICladeToPHatPLV(NNIClade clade_type) {
  switch (clade_type) {
    case NNIClade::ParentSister:
      return KeyIndex::Parent_PHatSister;
    case NNIClade::ChildLeft:
      return KeyIndex::Child_PHatLeft;
    case NNIClade::ChildRight:
      return KeyIndex::Child_PHatRight;
    default:
      Failwith("Given NNIClade has no associated KeyIndex.");
  }
};

NNIEngine::KeyIndexPairArray NNIEngine::BuildKeyIndexTypePairsFromPreNNIToPostNNI(
    const NNIOperation &pre_nni, const NNIOperation &post_nni) {
  // Find mapping from clades in pre-NNI to NNI.
  const auto nni_clade_map =
      NNIOperation::BuildNNICladeMapFromPreNNIToNNI(pre_nni, post_nni);
  KeyIndexPairArray key_idx_pair_array;
  size_t nni_clade_count = 0;
  for (const auto pre_nni_clade_type :
       {NNIClade::ParentSister, NNIClade::ChildLeft, NNIClade::ChildRight}) {
    const auto post_nni_clade_type = nni_clade_map[pre_nni_clade_type];
    key_idx_pair_array[nni_clade_count] = {NNICladeToPHatPLV(pre_nni_clade_type),
                                           NNICladeToPHatPLV(post_nni_clade_type)};
    nni_clade_count++;
  }
  return key_idx_pair_array;
}

NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForNNI(
    const NNIOperation &nni, const size_t node_count) const {
  return NNIEngine::BuildKeyIndexMapForNNI(nni, GetGraftDAG(), node_count);
}

NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForPostNNIViaReferencePreNNI(
    const NNIOperation &pre_nni, const NNIOperation &post_nni,
    const NNIEngine::KeyIndexMap &pre_key_idx) {
  return NNIEngine::BuildKeyIndexMapForPostNNIViaReferencePreNNI(
      pre_nni, post_nni, pre_key_idx, GetGraftDAG());
}

template <typename DAGType>
NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForNNI(const NNIOperation &nni,
                                                         const DAGType &dag,
                                                         const size_t node_count) {
  // Find NNI nodes.
  const auto parent_id = dag.GetDAGNodeId(nni.GetParent());
  const auto child_id = dag.GetDAGNodeId(nni.GetChild());
  const bool is_left_clade_sister = nni.WhichParentCladeIsFocalClade();
  // Find key indices for NNI.
  KeyIndexMap key_idx_map;
  key_idx_map.fill(NoId);
  key_idx_map[KeyIndex::Parent_Id] = parent_id.value_;
  key_idx_map[KeyIndex::Child_Id] = child_id.value_;
  key_idx_map[KeyIndex::Edge] = dag.GetEdgeIdx(parent_id, child_id).value_;
  key_idx_map[KeyIndex::Parent_RHat] =
      PLVHandler::GetPVIndex(PLVType::RHat, parent_id, node_count);
  key_idx_map[KeyIndex::Parent_RFocal] = PLVHandler::GetPVIndex(
      PLVHandler::RPLVType(!is_left_clade_sister), parent_id, node_count);
  key_idx_map[KeyIndex::Parent_PHatSister] = PLVHandler::GetPVIndex(
      PLVHandler::PPLVType(is_left_clade_sister), parent_id, node_count);
  key_idx_map[KeyIndex::Child_P] =
      PLVHandler::GetPVIndex(PLVType::P, child_id, node_count);
  key_idx_map[KeyIndex::Child_PHatLeft] =
      PLVHandler::GetPVIndex(PLVType::PHatLeft, child_id, node_count);
  key_idx_map[KeyIndex::Child_PHatRight] =
      PLVHandler::GetPVIndex(PLVType::PHatRight, child_id, node_count);

  return key_idx_map;
}
// Explicit Instantiation
template NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForNNI(
    const NNIOperation &nni, const GPDAG &dag, const size_t node_count);
template NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForNNI(
    const NNIOperation &nni, const GraftDAG &dag, const size_t node_count);

template <typename DAGType>
NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForPostNNIViaReferencePreNNI(
    const NNIOperation &pre_nni, const NNIOperation &post_nni,
    const NNIEngine::KeyIndexMap &pre_key_idx, const DAGType &dag) {
  // Unpopulated key indices will left as NoId.
  NodeId parent_id = dag.GetDAGNodeId(post_nni.GetParent());
  NodeId child_id = dag.GetDAGNodeId(post_nni.GetChild());
  KeyIndexMap post_key_idx;
  post_key_idx.fill(NoId);
  post_key_idx[KeyIndex::Parent_Id] = parent_id.value_;
  post_key_idx[KeyIndex::Child_Id] = child_id.value_;
  post_key_idx[KeyIndex::Edge] = dag.GetEdgeIdx(parent_id, child_id).value_;

  // Array for mapping from pre-NNI plvs to post-NNI plvs.
  const auto key_map = BuildKeyIndexTypePairsFromPreNNIToPostNNI(pre_nni, post_nni);

  // Set NNI plvs to their corresponding Pre-NNI plvs.
  post_key_idx[KeyIndex::Parent_RHat] = pre_key_idx[KeyIndex::Parent_RHat];
  for (const auto &[pre_key_type, post_key_type] : key_map) {
    post_key_idx[post_key_type] = pre_key_idx[pre_key_type];
  }

  return post_key_idx;
}
// Explicit Instantiation
template NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForPostNNIViaReferencePreNNI(
    const NNIOperation &pre_nni, const NNIOperation &post_nni,
    const NNIEngine::KeyIndexMap &pre_key_idx, const GPDAG &dag);
template NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForPostNNIViaReferencePreNNI(
    const NNIOperation &pre_nni, const NNIOperation &post_nni,
    const NNIEngine::KeyIndexMap &pre_key_idx, const GraftDAG &dag);

// ** Scoring via NNI Likelihood

NNIEngine::KeyIndexMapPair NNIEngine::PassGPEngineDataFromPreNNIToPostNNIViaCopy(
    const NNIOperation &pre_nni, const NNIOperation &post_nni) {
  // Find data in pre-NNI.
  const auto pre_key_idx =
      BuildKeyIndexMapForNNI(pre_nni, GetGraftDAG().NodeCount() - 1);
  // Find data in NNI.
  const auto &post_key_idx =
      BuildKeyIndexMapForNNI(post_nni, GetGraftDAG().NodeCount() - 1);

  // Array for mapping from pre-NNI plvs to post-NNI plvs.
  const auto key_type_map =
      BuildKeyIndexTypePairsFromPreNNIToPostNNI(pre_nni, post_nni);
  // Copy over pre-NNI plvs to NNI plvs.
  GetGPEngine().CopyPLVData(pre_key_idx[KeyIndex::Parent_RHat],
                            post_key_idx[KeyIndex::Parent_RHat]);
  for (const auto &[pre_key_type, post_key_type] : key_type_map) {
    GetGPEngine().CopyPLVData(pre_key_idx[post_key_type], post_key_idx[pre_key_type]);
  }

  // Copy over associated node data.
  for (const auto id_type : {KeyIndex::Parent_Id, KeyIndex::Child_Id}) {
    const auto pre_node_id = NodeId(pre_key_idx[id_type]);
    const auto post_node_id = NodeId(post_key_idx[id_type]);
    GetGPEngine().CopyNodeData(pre_node_id, post_node_id);
  }
  // Copy over central edge data.
  for (const auto idx_type : {KeyIndex::Edge}) {
    const auto pre_edge_idx = EdgeId(pre_key_idx[idx_type]);
    const auto post_edge_idx = EdgeId(post_key_idx[idx_type]);
    GetGPEngine().CopyGPCSPData(pre_edge_idx, post_edge_idx);
  }

  // Copy over associated non-central edge data.
  // Gather common ancestors and descendents of Pre-NNI and Post-NNI.
  for (const auto key_index : {KeyIndex::Parent_Id, KeyIndex::Child_Id}) {
    const auto pre_node = GetGraftDAG().GetDAGNode(NodeId(pre_key_idx[key_index]));
    const auto post_node = GetGraftDAG().GetDAGNode(NodeId(post_key_idx[key_index]));
    for (const auto direction : {Direction::Rootward, Direction::Leafward}) {
      // Ignore parents of child node.
      if ((key_index == KeyIndex::Child_Id) && (direction == Direction::Rootward)) {
        continue;
      }
      for (const auto is_focal_clade : {true, false}) {
        // Ignore focal children of parent node.
        if ((key_index == KeyIndex::Parent_Id) && (direction == Direction::Leafward) &&
            is_focal_clade) {
          continue;
        }
        SubsplitClade prenni_clade = (is_focal_clade ? pre_nni.WhichCladeIsFocal()
                                                     : pre_nni.WhichCladeIsSister());
        for (const auto &adj_node_id : pre_node.GetNeighbors(direction, prenni_clade)) {
          // If edge from Pre-NNI also exists in Post-NNI, copy data over.
          if (GetGraftDAG().ContainsEdge(post_node.Id(), NodeId(adj_node_id))) {
            const auto pre_edge_idx =
                GetGraftDAG().GetEdgeIdx(pre_node.Id(), NodeId(adj_node_id));
            const auto post_edge_idx =
                GetGraftDAG().GetEdgeIdx(post_node.Id(), NodeId(adj_node_id));
            GetGPEngine().CopyGPCSPData(pre_edge_idx, post_edge_idx);
          }
        }
      }
    }
  }

  return {pre_key_idx, post_key_idx};
}

NNIEngine::KeyIndexMapPair NNIEngine::PassGPEngineDataFromPreNNIToPostNNIViaReference(
    const NNIOperation &pre_nni, const NNIOperation &post_nni, const size_t nni_count,
    const bool use_unique_temps) {
  // Find data in pre-NNI.
  KeyIndexMap pre_key_idx =
      BuildKeyIndexMapForNNI(pre_nni, GetGPDAG().NodeCountWithoutDAGRoot());
  // Find data in post-NNI.
  KeyIndexMap post_key_idx =
      BuildKeyIndexMapForPostNNIViaReferencePreNNI(pre_nni, post_nni, pre_key_idx);

  // Assign temporary place in engine to store intermediate values.
  size_t temp_offset_1, temp_offset_2;
  if (use_unique_temps) {
    temp_offset_1 = nni_count * 2;
    temp_offset_2 = (nni_count * 2) + 1;
  } else {
    temp_offset_1 = 0;
    temp_offset_2 = 1;
  }
  post_key_idx[KeyIndex::Parent_RFocal] = GetGPEngine().GetSparePLVIndex(temp_offset_1);
  post_key_idx[KeyIndex::Child_P] = GetGPEngine().GetSparePLVIndex(temp_offset_2);
  post_key_idx[KeyIndex::Edge] = GetGPEngine().GetSpareGPCSPIndex(nni_count);

  // Copy over central edge data.
  for (const auto idx_type : {KeyIndex::Edge}) {
    const auto pre_edge_idx = EdgeId(pre_key_idx[idx_type]);
    const auto post_edge_idx = EdgeId(post_key_idx[idx_type]);
    GetGPEngine().CopyGPCSPData(pre_edge_idx, post_edge_idx);
  }

  return {pre_key_idx, post_key_idx};
}

GPOperationVector NNIEngine::BuildGPOperationsForNNILikelihood(
    const NNIOperation &nni, const KeyIndexMap &nni_key_idx) const {
  GPOperationVector ops;
  // p(parent_focal) = rhat(parent) \circ phat(parent_sister)
  ops.push_back(GPOperations::Multiply{nni_key_idx[KeyIndex::Parent_RFocal],
                                       nni_key_idx[KeyIndex::Parent_RHat],
                                       nni_key_idx[KeyIndex::Parent_PHatSister]});
  // p(child) = phat(child_left) \circ phat(child_right)
  ops.push_back(GPOperations::Multiply{nni_key_idx[KeyIndex::Child_P],
                                       nni_key_idx[KeyIndex::Child_PHatLeft],
                                       nni_key_idx[KeyIndex::Child_PHatRight]});
  // compute likelihood on edge between r(parent_focal) and p(child)
  ops.push_back(GPOperations::Likelihood{nni_key_idx[KeyIndex::Edge],
                                         nni_key_idx[KeyIndex::Parent_RFocal],
                                         nni_key_idx[KeyIndex::Child_P]});
  return ops;
}

GPOperationVector NNIEngine::BuildGPOperationsForAdjacentNNILikelihoods(
    const bool via_reference) {
  GPOperationVector operations = GPOperationVector();
  GrowGPEngineForAdjacentNNILikelihoods(via_reference);
  size_t nni_count = 0;
  for (const auto &nni : GetAdjacentNNIs()) {
    const auto pre_nni = FindNNINeighborInDAG(nni);
    const auto [prenni_key_idx, nni_key_idx] =
        (via_reference ? PassGPEngineDataFromPreNNIToPostNNIViaReference(
                             pre_nni, nni, nni_count, false)
                       : PassGPEngineDataFromPreNNIToPostNNIViaCopy(pre_nni, nni));
    std::ignore = prenni_key_idx;
    GPOperations::AppendGPOperations(
        operations, BuildGPOperationsForNNILikelihood(nni, nni_key_idx));
    nni_count++;
  }
  return operations;
}

void NNIEngine::GrowGPEngineForAdjacentNNILikelihoods(const bool via_reference,
                                                      const bool use_unique_temps) {
  if (via_reference) {
    if (use_unique_temps) {
      GetGPEngine().GrowSparePLVs(2 * GetAdjacentNNICount());
    } else {
      GetGPEngine().GrowSparePLVs(2);
    }
    GetGPEngine().GrowSpareGPCSPs(GetAdjacentNNICount());
  } else {
    GetGPEngine().GrowPLVs(GetGraftDAG().NodeCountWithoutDAGRoot());
    GetGPEngine().GrowGPCSPs(GetGraftDAG().EdgeCountWithLeafSubsplits());
  }
}

void NNIEngine::ScoreAdjacentNNIsByGPLikelihood() {
  // Compute Likelihoods
  const GPOperationVector operations = BuildGPOperationsForAdjacentNNILikelihoods(true);
  GetGPEngine().ProcessOperations(operations);
  // Retrieve results from GPEngine and store in Scored NNIs.
  const auto proposed_nni_likelihoods =
      GetGPEngine().GetSparePerGPCSPLogLikelihoods(0, GetAdjacentNNICount());
  size_t nni_count = 0;
  for (const auto &nni : GetAdjacentNNIs()) {
    scored_nnis_[nni] = proposed_nni_likelihoods[nni_count];
    nni_count++;
  }
}

// ** Scoring via TP Likelihood

void NNIEngine::ScoreAdjacentNNIsByTPLikelihood() {
  // Retrieve results from TPEngine and store in Scored NNIs.
  for (const auto &nni : GetAdjacentNNIs()) {
    const auto pre_nni = FindNNINeighborInDAG(nni);
    scored_nnis_[nni] = GetTPEngine().GetTopTreeLikelihoodWithProposedNNI(nni, pre_nni);
  }
}

// ** Scoring via TP Parsimony

void NNIEngine::ScoreAdjacentNNIsByTPParsimony() {
  // !ISSUE: Implement Parsimony Score
  Failwith("Error: No implementation.");
}

// ** DAG & Engine Maintenance

void NNIEngine::InitGPEngine() {
  if (gp_engine_) {
    GetGPEngine().GrowPLVs(dag_.NodeCountWithoutDAGRoot());
    GetGPEngine().GrowGPCSPs(dag_.EdgeCountWithLeafSubsplits());
    GetGPEngine().ProcessOperations(dag_.PopulatePLVs());
    GetGPEngine().ProcessOperations(dag_.ComputeLikelihoods());
  }
  if (tp_engine_) {
    GetTPEngine().GrowNodeData(dag_.NodeCount());
    GetTPEngine().GrowEdgeData(dag_.EdgeCountWithLeafSubsplits());
  }
}

void NNIEngine::PrepGPEngineForLikelihoods() {
  GetGPEngine().ProcessOperations(dag_.PopulatePLVs());
  GetGPEngine().ProcessOperations(dag_.ComputeLikelihoods());
}

void NNIEngine::AddAcceptedNNIsToDAG() {
  // Make sure all grafted nodes have been removed from the DAG
  graft_dag_->RemoveAllGrafts();
  // Initialize reindexers for remapping after adding nodes.
  node_reindexer_ = Reindexer::IdentityReindexer(dag_.NodeCount());
  edge_reindexer_ = Reindexer::IdentityReindexer(dag_.EdgeCountWithLeafSubsplits());
  for (const auto &nni : accepted_nnis_) {
    auto mods = dag_.AddNodePair(nni);
    node_reindexer_.ComposeWith(mods.node_reindexer);
    edge_reindexer_.ComposeWith(mods.edge_reindexer);
  }
  // Remove DAGRoot from node reindexing.
  node_reindexer_ = node_reindexer_.RemoveNewIndex(dag_.GetDAGRootNodeId().value_);
  // Grow GPEngine to fit accepted NNIs.
  Assert(dag_.NodeCountWithoutDAGRoot() == node_reindexer_.size(),
         "Node reindexer is the wrong size.");
  Assert(dag_.EdgeCountWithLeafSubsplits() == edge_reindexer_.size(),
         "Edge reindexer is the wrong size.");
  // Resize engines.
  if (gp_engine_) {
    GetGPEngine().GrowPLVs(dag_.NodeCountWithoutDAGRoot(), node_reindexer_);
    GetGPEngine().GrowGPCSPs(dag_.EdgeCountWithLeafSubsplits(), edge_reindexer_);
  }
  if (tp_engine_) {
    GetTPEngine().GrowNodeData(dag_.NodeCount(), node_reindexer_);
    GetTPEngine().GrowEdgeData(dag_.EdgeCountWithLeafSubsplits(), edge_reindexer_);
  }
}

void NNIEngine::GraftAdjacentNNIsToDAG() {
  for (const auto &nni : GetAdjacentNNIs()) {
    graft_dag_->AddNodePair(nni);
  }
}

void NNIEngine::RemoveAllGraftedNNIsFromDAG() { graft_dag_->RemoveAllGrafts(); }

// ** NNI Maintenance

void NNIEngine::SyncAdjacentNNIsWithDAG() {
  adjacent_nnis_.clear();
  // Only real node pairs are viable NNIs.
  dag_.IterateOverRealNodes([this](SubsplitDAGNode node) {
    dag_.IterateOverParentAndChildAndLeafwardEdges(
        node, [this](const NodeId parent_id, const bool is_edge_on_left,
                     const NodeId child_id, const EdgeId edge_idx) {
          // Only internal node pairs are viable NNIs.
          const Bitset &parent_bitset = dag_.GetDAGNode(parent_id).GetBitset();
          const Bitset &child_bitset = dag_.GetDAGNode(child_id).GetBitset();
          if (!(parent_bitset.SubsplitIsUCA() || child_bitset.SubsplitIsLeaf())) {
            SafeAddOutputNNIsToAdjacentNNIs(parent_bitset, child_bitset,
                                            is_edge_on_left);
          }
        });
  });
}

void NNIEngine::UpdateAdjacentNNIsAfterDAGAddNodePair(const NNIOperation &nni) {
  UpdateAdjacentNNIsAfterDAGAddNodePair(nni.parent_, nni.child_);
}

void NNIEngine::UpdateAdjacentNNIsAfterDAGAddNodePair(const Bitset &parent_bitset,
                                                      const Bitset &child_bitset) {
  const auto parent_id = dag_.GetDAGNodeId(parent_bitset);
  const auto child_id = dag_.GetDAGNodeId(child_bitset);
  // Every new edge added is a potential new NNI.
  // Iterate over the parent and child node of the new pair.
  for (const auto &node_id : {parent_id, child_id}) {
    // Get nodes adjacent to current node from both left and right edges.
    for (const bool is_edge_leafward : {true, false}) {
      // Get nodes adjacent to current node from both leafward and rootward
      // directions.
      for (const bool is_edge_on_left : {true, false}) {
        auto adjacent_node_ids = dag_.GetDAGNode(node_id).GetLeafwardOrRootward(
            is_edge_leafward, is_edge_on_left);
        AddAllNNIsFromNodeVectorToAdjacentNNIs(node_id, adjacent_node_ids,
                                               is_edge_on_left, is_edge_leafward);
      }
    }
  }
  // Remove the pair that was just added to the DAG from NNI Set.
  NNIOperation new_nni = NNIOperation(parent_bitset, child_bitset);
  adjacent_nnis_.erase(new_nni);
}

void NNIEngine::AddAllNNIsFromNodeVectorToAdjacentNNIs(
    const NodeId node_id, const SizeVector &adjacent_node_ids,
    const bool is_edge_on_left, const bool is_edge_leafward) {
  Bitset node_bitset = dag_.GetDAGNode(node_id).GetBitset();
  // Determine whether node_id corresponds to parent or child of the pair.
  // Add every edge's NNI to NNI Set.
  // If edges are leafward, node_id is the parent to all vector nodes.
  // If edges are rootward, node_id is the child to all vector nodes.
  if (is_edge_leafward) {
    const Bitset &parent_bitset = node_bitset;
    for (const auto &adjacent_node_id : adjacent_node_ids) {
      const Bitset child_bitset = dag_.GetDAGNode(NodeId(adjacent_node_id)).GetBitset();
      SafeAddOutputNNIsToAdjacentNNIs(parent_bitset, child_bitset, is_edge_on_left);
    }
  } else {
    const Bitset &child_bitset = node_bitset;
    for (const auto &adjacent_node_id : adjacent_node_ids) {
      const Bitset parent_bitset =
          dag_.GetDAGNode(NodeId(adjacent_node_id)).GetBitset();
      SafeAddOutputNNIsToAdjacentNNIs(parent_bitset, child_bitset, is_edge_on_left);
    }
  }
}

void NNIEngine::SafeAddOutputNNIsToAdjacentNNIs(const Bitset &parent_bitset,
                                                const Bitset &child_bitset,
                                                const bool is_edge_on_left) {
  // Soft assert that parent is not the root and child is not a leaf.
  if (parent_bitset.SubsplitIsUCA() || child_bitset.SubsplitIsLeaf()) {
    return;
  }
  // Input pair is in the DAG, so remove it from the Set if it exists.
  adjacent_nnis_.erase({parent_bitset, child_bitset});
  // Add NNI for right clade swap and left clade swap.
  for (bool is_swap_with_left_child : {true, false}) {
    bool is_in_dag = false;
    const auto new_nni = NNIOperation::NNIOperationFromNeighboringSubsplits(
        parent_bitset, child_bitset, is_swap_with_left_child, !is_edge_on_left);
    // If DAG already contains output parent and child nodes, and an edge between
    // them, then don't add it to the adjacent_nnis.
    if (dag_.ContainsNode(new_nni.parent_) && dag_.ContainsNode(new_nni.child_)) {
      const auto parent_id = dag_.GetDAGNodeId(new_nni.parent_);
      const auto child_id = dag_.GetDAGNodeId(new_nni.child_);
      is_in_dag = dag_.ContainsEdge(parent_id, child_id);
    }
    if (!is_in_dag) {
      adjacent_nnis_.insert(new_nni);
    }
  }
}

void NNIEngine::AddScoreForNNI(const NNIOperation &nni, const double score) {
  SafeInsert(scored_nnis_, nni, score);
};

double NNIEngine::GetScoreForNNI(const NNIOperation &nni) const {
  Assert(scored_nnis_.find(nni) != scored_nnis_.end(),
         "NNIEngine::GetScoreForNNI(): Score does not exist for requested NNI.");
  const double score = scored_nnis_.at(nni);
  return score;
}

void NNIEngine::UpdateAdjacentNNIs(const bool reevaluate_rejected_nnis) {
  adjacent_nnis_.clear();
  if (reevaluate_rejected_nnis) {
    adjacent_nnis_.insert(rejected_nnis_.begin(), rejected_nnis_.end());
  }
  for (const auto &nni : GetAcceptedNNIs()) {
    const auto focal_clade =
        Bitset::SubsplitIsChildOfWhichParentClade(nni.parent_, nni.child_);
    const bool is_edge_on_left = (focal_clade == SubsplitClade::Left);
    SafeAddOutputNNIsToAdjacentNNIs(nni.parent_, nni.child_, is_edge_on_left);
  }
}

void NNIEngine::UpdateAcceptedNNIs(const bool save_past_nnis) {
  if (save_past_nnis) {
    accepted_past_nnis_.insert(accepted_nnis_.begin(), accepted_nnis_.end());
  }
  accepted_nnis_.clear();
}

void NNIEngine::UpdateRejectedNNIs(const bool save_past_nnis) {
  if (save_past_nnis) {
    rejected_past_nnis_.insert(rejected_nnis_.begin(), rejected_nnis_.end());
  }
  rejected_nnis_.clear();
}

void NNIEngine::UpdateScoredNNIs(const bool save_past_nnis) {
  if (save_past_nnis) {
    scored_past_nnis_.insert(scored_nnis_.begin(), scored_nnis_.end());
  }
  scored_nnis_.clear();
}

void NNIEngine::ResetAllNNIs() {
  adjacent_nnis_.clear();
  accepted_nnis_.clear();
  accepted_past_nnis_.clear();
  rejected_nnis_.clear();
  rejected_past_nnis_.clear();
}

// ** Miscellaneous

NNIOperation NNIEngine::FindNNINeighborInDAG(const NNIOperation &nni) {
  for (const auto which_swap : {true, false}) {
    const auto &swapped_nni = nni.NNIOperationFromNeighboringSubsplits(which_swap);
    if (dag_.ContainsNode(swapped_nni.parent_) &&
        dag_.ContainsNode(swapped_nni.child_)) {
      auto parent_id = dag_.GetDAGNodeId(swapped_nni.parent_);
      auto child_id = dag_.GetDAGNodeId(swapped_nni.child_);
      if (parent_id > dag_.NodeCount() || child_id > dag_.NodeCount()) {
        std::cout << "ERROR: found NNI Neighbor from the GraftDAG!" << std::endl;
      }
      return swapped_nni;
    }
  }
  Failwith("NNIOperation has no neighbors found in the DAG.");
}
