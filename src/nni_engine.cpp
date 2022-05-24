// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "nni_engine.hpp"

#include "gp_engine.hpp"
#include "gp_dag.hpp"

#include "bitset.hpp"
#include "subsplit_dag.hpp"
#include "nni_operation.hpp"
#include "graft_dag.hpp"
#include "sugar.hpp"

using PLVType = PLVHandler::PLVType;

NNIEngine::NNIEngine(GPDAG &dag, GPEngine &gp_engine)
    : dag_(dag), graft_dag_(std::make_unique<GraftDAG>(dag)), gp_engine_(gp_engine) {}

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
  FilterInit();
}

void NNIEngine::RunMainLoop() {
  // (1) Add all adjacent NNIs to the GraftDAG.
  GPEngineUpdate();
  GraftAdjacentNNIsToDAG();
  // (2) Select whether to accept or reject adjacent NNIs via filter.
  FilterPreUpdate();
  FilterProcessAdjacentNNIs();
  FilterPostUpdate();

  sweep_count_++;
}

void NNIEngine::RunPostLoop() {
  // (3) Add accepted NNIs to permanent DAG.
  AddAcceptedNNIsToDAG();
  // (4) Update Adjacent NNIs to reflect added NNI.
  UpdateAdjacentNNIs();
  // (5) Reset Accepted NNIs and GraftDAG and save results.
  UpdateAcceptedNNIs();
  UpdateRejectedNNIs();
}

// ** Static Filter Functions

void NNIEngine::NoInit(NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                       GraftDAG &this_graft_dag) {}

void NNIEngine::NoUpdate(NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                         GraftDAG &this_graft_dag) {}

double NNIEngine::GetBestPreNNILikelihood(const NNIOperation &nni) {
  size_t prenni_edge_idx = NoId;
  double prenni_likelihood = -INFINITY;
  SizeVector prenni_idxs;
  DoubleVector prenni_likelihoods;
  for (const SubsplitClade clade_to_swap :
       {SubsplitClade::Right, SubsplitClade::Left}) {
    NNIOperation pre_nni = nni.NNIOperationFromNeighboringSubsplits(
        (clade_to_swap == SubsplitClade::Right));
    if (dag_.ContainsNode(pre_nni.GetParent()) &&
        dag_.ContainsNode(pre_nni.GetChild())) {
      const auto &parent_id = dag_.GetDAGNodeId(pre_nni.GetParent());
      const auto &child_id = dag_.GetDAGNodeId(pre_nni.GetChild());
      if (dag_.ContainsEdge(parent_id, child_id)) {
        const auto &edge_idx = dag_.GetEdgeIdx(parent_id, child_id);
        const auto likelihood = GetDAGEdgeLikelihood(edge_idx);
        if (prenni_likelihood < likelihood) {
          prenni_likelihood = likelihood;
          prenni_edge_idx = edge_idx;
        }
        prenni_idxs.push_back(edge_idx);
        prenni_likelihoods.push_back(likelihood);
      }
    }
  }
  Assert(prenni_edge_idx != NoId,
         "Pre-NNI could not be found in DAG based on Post-NNI.");
  return prenni_likelihood;
}

double NNIEngine::GetPostNNILikelihood(const NNIOperation &nni) {
  return scored_nnis_[nni];
}

void NNIEngine::UpdateScoreAdjacentNNIsByLikelihood(NNIEngine &this_nni_engine,
                                                    GPEngine &this_gp_engine,
                                                    GraftDAG &this_graft_dag) {
  this_nni_engine.ScoreAdjacentNNIsByLikelihood();
}

bool NNIEngine::FilterPostNNIBetterThanPreNNI(NNIEngine &this_nni_engine,
                                              GPEngine &this_gp_engine,
                                              GraftDAG &this_graft_dag,
                                              const NNIOperation &nni) {
  const double postnni_likelihood = this_nni_engine.GetPostNNILikelihood(nni);
  const double prenni_likelihood = this_nni_engine.GetBestPreNNILikelihood(nni);
  const double diff = postnni_likelihood - prenni_likelihood;
  const bool accept = (diff > 0.0);
  return accept;
}

// ** NNI Likelihoods

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

NNIEngine::KeyIndexMapPair NNIEngine::PassDataFromPreNNIToPostNNIViaCopy(
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
  gp_engine_.CopyPLVData(pre_key_idx[KeyIndex::Parent_RHat],
                         post_key_idx[KeyIndex::Parent_RHat]);
  for (const auto &[pre_key_type, post_key_type] : key_type_map) {
    gp_engine_.CopyPLVData(pre_key_idx[post_key_type], post_key_idx[pre_key_type]);
  }

  // Copy over associated node data.
  for (const auto id_type : {KeyIndex::Parent_Id, KeyIndex::Child_Id}) {
    const auto pre_node_id = pre_key_idx[id_type];
    const auto post_node_id = post_key_idx[id_type];
    gp_engine_.CopyNodeData(pre_node_id, post_node_id);
  }
  // Copy over central edge data.
  for (const auto idx_type : {KeyIndex::Edge}) {
    const auto pre_edge_idx = pre_key_idx[idx_type];
    const auto post_edge_idx = post_key_idx[idx_type];
    gp_engine_.CopyGPCSPData(pre_edge_idx, post_edge_idx);
  }

  // Copy over associated non-central edge data.
  // Gather common ancestors and descendents of Pre-NNI and Post-NNI.
  for (const auto key_index : {KeyIndex::Parent_Id, KeyIndex::Child_Id}) {
    const auto pre_node = GetGraftDAG().GetDAGNode(pre_key_idx[key_index]);
    const auto post_node = GetGraftDAG().GetDAGNode(post_key_idx[key_index]);
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
          if (GetGraftDAG().ContainsEdge(post_node.Id(), adj_node_id)) {
            const auto pre_edge_idx =
                GetGraftDAG().GetEdgeIdx(pre_node.Id(), adj_node_id);
            const auto post_edge_idx =
                GetGraftDAG().GetEdgeIdx(post_node.Id(), adj_node_id);
            gp_engine_.CopyGPCSPData(pre_edge_idx, post_edge_idx);
          }
        }
      }
    }
  }
  return {pre_key_idx, post_key_idx};
}

NNIEngine::KeyIndexMapPair NNIEngine::PassDataFromPreNNIToPostNNIViaReference(
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
  post_key_idx[KeyIndex::Parent_RFocal] = gp_engine_.GetTempPLVIndex(temp_offset_1);
  post_key_idx[KeyIndex::Child_P] = gp_engine_.GetTempPLVIndex(temp_offset_2);
  post_key_idx[KeyIndex::Edge] = gp_engine_.GetTempGPCSPIndex(nni_count);

  // Copy over central edge data.
  for (const auto idx_type : {KeyIndex::Edge}) {
    const auto pre_edge_idx = pre_key_idx[idx_type];
    const auto post_edge_idx = post_key_idx[idx_type];
    gp_engine_.CopyGPCSPData(pre_edge_idx, post_edge_idx);
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
  GrowEngineForAdjacentNNILikelihoods(via_reference);
  size_t nni_count = 0;
  for (const auto &nni : GetAdjacentNNIs()) {
    const auto pre_nni = FindNNINeighborInDAG(nni);
    const auto [prenni_key_idx, nni_key_idx] =
        (via_reference
             ? PassDataFromPreNNIToPostNNIViaReference(pre_nni, nni, nni_count, false)
             : PassDataFromPreNNIToPostNNIViaCopy(pre_nni, nni));
    std::ignore = prenni_key_idx;
    GPOperations::AppendGPOperations(
        operations, BuildGPOperationsForNNILikelihood(nni, nni_key_idx));
    nni_count++;
  }
  return operations;
}

void NNIEngine::GrowEngineForAdjacentNNILikelihoods(const bool via_reference,
                                                    const bool use_unique_temps) {
  if (via_reference) {
    if (use_unique_temps) {
      gp_engine_.GrowTempPLVs(2 * GetAdjacentNNICount());
    } else {
      gp_engine_.GrowTempPLVs(2);
    }
    gp_engine_.GrowTempGPCSPs(GetAdjacentNNICount());
  } else {
    gp_engine_.GrowPLVs(GetGraftDAG().NodeCountWithoutDAGRoot());
    gp_engine_.GrowGPCSPs(GetGraftDAG().EdgeCountWithLeafSubsplits());
  }
}

void NNIEngine::ScoreAdjacentNNIsByLikelihood() {
  // Compute Likelihoods
  const GPOperationVector operations = BuildGPOperationsForAdjacentNNILikelihoods(true);
  gp_engine_.ProcessOperations(operations);
  // Retrieve results from GPEngine and store in Scored NNIs.
  size_t nni_count = 0;
  const auto likelihoods =
      gp_engine_.GetTempPerGPCSPLogLikelihoods(0, GetAdjacentNNICount());
  for (const auto &nni : GetAdjacentNNIs()) {
    scored_nnis_[nni] = likelihoods[nni_count];
    nni_count++;
  }
}

// ** Runner Subroutines

void NNIEngine::FilterInit() { filter_init_fn_(*this, gp_engine_, GetGraftDAG()); }

void NNIEngine::FilterPreUpdate() {
  filter_pre_update_fn_(*this, gp_engine_, GetGraftDAG());
}

void NNIEngine::FilterProcessAdjacentNNIs() {
  for (const auto &nni : GetAdjacentNNIs()) {
    const bool accept_nni = filter_process_fn_(*this, gp_engine_, GetGraftDAG(), nni);
    if (accept_nni) {
      accepted_nnis_.insert(nni);
    } else {
      rejected_nnis_.insert(nni);
    }
  }
}

void NNIEngine::FilterPostUpdate() {
  filter_post_update_fn_(*this, gp_engine_, GetGraftDAG());
}

// ** Set Filtering Scheme

void NNIEngine::SetFilteringCustom(StaticFilterInitFunction init_fn,
                                   StaticFilterUpdateFunction pre_update_fn,
                                   StaticFilterProcessFunction process_fn,
                                   StaticFilterUpdateFunction post_update_fn) {
  filter_init_fn_ = init_fn;
  filter_pre_update_fn_ = pre_update_fn;
  filter_process_fn_ = process_fn;
  filter_post_update_fn_ = post_update_fn;
}

void NNIEngine::SetFilteringNone() {
  filter_init_fn_ = NoInit;
  filter_pre_update_fn_ = NoUpdate;
  filter_process_fn_ = NoFilter<true>;
  filter_post_update_fn_ = NoUpdate;
}

void NNIEngine::SetFilteringHillclimb() {
  filter_init_fn_ = NoInit;
  filter_pre_update_fn_ = UpdateScoreAdjacentNNIsByLikelihood;
  filter_process_fn_ = FilterPostNNIBetterThanPreNNI;
  filter_post_update_fn_ = NoUpdate;
}

void NNIEngine::SetFilteringLossThreshold(const double loss_tolerance) {
  filter_init_fn_ = NoInit;
  filter_pre_update_fn_ = UpdateScoreAdjacentNNIsByLikelihood;
  filter_process_fn_ = [loss_tolerance](
                           NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                           GraftDAG &this_graft_dag, const NNIOperation &nni) -> bool {
    const double postnni_likelihood = this_nni_engine.GetPostNNILikelihood(nni);
    const double prenni_likelihood = this_nni_engine.GetBestPreNNILikelihood(nni);
    const double diff = postnni_likelihood - prenni_likelihood;
    return diff >= loss_tolerance;
  };
  filter_post_update_fn_ = NoUpdate;
}

void NNIEngine::SetFilteringCutoffThreshold(const double cutoff_threshold) {
  filter_init_fn_ = NoInit;
  filter_pre_update_fn_ = UpdateScoreAdjacentNNIsByLikelihood;
  filter_process_fn_ = [cutoff_threshold](
                           NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                           GraftDAG &this_graft_dag, const NNIOperation &nni) -> bool {
    const double postnni_likelihood = this_nni_engine.GetPostNNILikelihood(nni);
    return postnni_likelihood >= cutoff_threshold;
  };
  filter_post_update_fn_ = NoUpdate;
}

// ** DAG & GPEngine Maintenance

void NNIEngine::GPEngineInit() { GPEngineUpdate(); }

void NNIEngine::GPEngineUpdate() {
  gp_engine_.GrowPLVs(dag_.NodeCountWithoutDAGRoot());
  gp_engine_.GrowGPCSPs(dag_.EdgeCountWithLeafSubsplits());
  gp_engine_.ProcessOperations(dag_.PopulatePLVs());
  if (do_optimitize_branch_lengths_) {
    gp_engine_.ProcessOperations(dag_.BranchLengthOptimization());
  }
  gp_engine_.ProcessOperations(dag_.ComputeLikelihoods());
  dag_likelihoods_ = std::make_unique<EigenVectorXd>(
      gp_engine_.GetPerGPCSPLogLikelihoods(0, dag_.EdgeCountWithLeafSubsplits()));
}

void NNIEngine::AddAcceptedNNIsToDAG() {
  // Make sure all grafted nodes have been removed from the DAG
  RemoveAllGraftedNNIsFromDAG();
  // Initialize reindexers for remapping after adding nodes.
  node_reindexer_ = Reindexer::IdentityReindexer(dag_.NodeCount());
  edge_reindexer_ = Reindexer::IdentityReindexer(dag_.EdgeCountWithLeafSubsplits());
  for (const auto &nni : accepted_nnis_) {
    auto mods = dag_.AddNodePair(nni);
    node_reindexer_.ComposeWith(mods.node_reindexer);
    edge_reindexer_.ComposeWith(mods.edge_reindexer);
  }
  // Remove DAGRoot from node reindexing.
  node_reindexer_ = node_reindexer_.RemoveNewIndex(dag_.GetDAGRootNodeId());
  // Grow GPEngine to fit accepted NNIs.
  Assert(dag_.NodeCountWithoutDAGRoot() == node_reindexer_.size(),
         "Node reindexer is the wrong size.");
  Assert(dag_.EdgeCountWithLeafSubsplits() == edge_reindexer_.size(),
         "Edge reindexer is the wrong size.");
  gp_engine_.GrowPLVs(dag_.NodeCountWithoutDAGRoot(), node_reindexer_);
  gp_engine_.GrowGPCSPs(dag_.EdgeCountWithLeafSubsplits(), edge_reindexer_);
}

void NNIEngine::GraftAdjacentNNIsToDAG() {
  RemoveAllGraftedNNIsFromDAG();
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
        node, [this](const size_t parent_id, const bool is_edge_on_left,
                     const size_t child_id, const size_t edge_idx) {
          // Only internal node pairs are viable NNIs.
          const Bitset &parent_bitset = dag_.GetDAGNode(parent_id).GetBitset();
          const Bitset &child_bitset = dag_.GetDAGNode(child_id).GetBitset();
          SafeAddOutputNNIsToAdjacentNNIs(parent_bitset, child_bitset, is_edge_on_left);
        });
  });
}

void NNIEngine::UpdateAdjacentNNIsAfterDAGAddNodePair(const NNIOperation &nni) {
  UpdateAdjacentNNIsAfterDAGAddNodePair(nni.parent_, nni.child_);
}

void NNIEngine::UpdateAdjacentNNIsAfterDAGAddNodePair(const Bitset &parent_bitset,
                                                      const Bitset &child_bitset) {
  const size_t parent_id = dag_.GetDAGNodeId(parent_bitset);
  const size_t child_id = dag_.GetDAGNodeId(child_bitset);
  // Every new edge added is a potential new NNI.
  // Iterate over the parent and child node of the new pair.
  for (const size_t &node_id : {parent_id, child_id}) {
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
    const size_t &node_id, const SizeVector &adjacent_node_ids,
    const bool is_edge_on_left, const bool is_edge_leafward) {
  Bitset node_bitset = dag_.GetDAGNode(node_id).GetBitset();
  // Determine whether node_id corresponds to parent or child of the pair.
  // Add every edge's NNI to NNI Set.
  // If edges are leafward, node_id is the parent to all vector nodes.
  // If edges are rootward, node_id is the child to all vector nodes.
  if (is_edge_leafward) {
    const Bitset &parent_bitset = node_bitset;
    for (const auto &adjacent_node_id : adjacent_node_ids) {
      const Bitset child_bitset = dag_.GetDAGNode(adjacent_node_id).GetBitset();
      SafeAddOutputNNIsToAdjacentNNIs(parent_bitset, child_bitset, is_edge_on_left);
    }
  } else {
    const Bitset &child_bitset = node_bitset;
    for (const auto &adjacent_node_id : adjacent_node_ids) {
      const Bitset parent_bitset = dag_.GetDAGNode(adjacent_node_id).GetBitset();
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
      const size_t parent_id = dag_.GetDAGNodeId(new_nni.parent_);
      const size_t child_id = dag_.GetDAGNodeId(new_nni.child_);
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

void NNIEngine::UpdateAdjacentNNIs() {
  adjacent_nnis_.clear();
  if (reevaluate_rejected_nnis_) {
    adjacent_nnis_.insert(rejected_nnis_.begin(), rejected_nnis_.end());
  }

  for (const auto &nni : GetAcceptedNNIs()) {
    for (const auto &node_subsplit : {nni.GetParent(), nni.GetChild()}) {
      const auto &node = dag_.GetDAGNode(dag_.GetDAGNodeId(node_subsplit));
      for (const auto &direction : {Direction::Rootward, Direction::Leafward}) {
        for (const auto &focal_clade : {SubsplitClade::Left, SubsplitClade::Right}) {
          const auto &neighbors = node.GetNeighbors(direction, focal_clade);
          for (const auto &adj_node_id : neighbors) {
            const auto &adj_node = dag_.GetDAGNode(adj_node_id);
            const auto &parent = (direction == Direction::Rootward) ? adj_node : node;
            const auto &child = (direction == Direction::Rootward) ? node : adj_node;
            const bool is_edge_on_left = (focal_clade == SubsplitClade::Left);
            SafeAddOutputNNIsToAdjacentNNIs(parent.GetBitset(), child.GetBitset(),
                                            is_edge_on_left);
          }
        }
      }
    }
  }
}

void NNIEngine::UpdateAcceptedNNIs() {
  if (keep_accepted_past_nnis_) {
    accepted_past_nnis_.insert(accepted_nnis_.begin(), accepted_nnis_.end());
  }
  accepted_nnis_.clear();
}

void NNIEngine::UpdateRejectedNNIs() {
  if (keep_rejected_past_nnis_) {
    rejected_past_nnis_.insert(rejected_nnis_.begin(), rejected_nnis_.end());
  }
  rejected_nnis_.clear();
}

void NNIEngine::UpdateScoredNNIs() {
  if (keep_scored_past_nnis_) {
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
      return swapped_nni;
    }
  }
  Failwith("NNIOperation has no neighbors found in the DAG.");
}
