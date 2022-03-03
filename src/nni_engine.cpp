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

NNIEngine::NNIEngine(GPDAG &dag, GPEngine &gp_engine)
    : dag_(dag),
      graft_dag_(std::make_unique<GraftDAG>(dag)),
      gp_engine_(gp_engine),
      run_count_(0) {}

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
  ClearAdjacentNNIs();
  SyncAdjacentNNIsWithDAG();
  // Initialize GPEngine for current DAG size
  ResizeGPEngineForDAG();
  PrepGPEngineForLikelihoods();
}

void NNIEngine::RunMainLoop() {
  // (1) Add all adjacent NNIs to the GraftDAG.
  GraftAdjacentNNIsToDAG();
  // (2) Compute marginal likelihoods for each adjacent NNI.
  EvaluateAdjacentNNIs();
  // (3) Select whether to accept or reject adjacent NNIs via filter.
  FilterAdjacentNNIsToAcceptedNNIs();
  // (4) Add accepted NNIs to permanent DAG.
  AddAcceptedNNIsToDAG();

  run_count_++;
}

void NNIEngine::RunPostLoop() {
  // (5) Update Adjacent NNIs to reflect added NNI.
  ClearAdjacentNNIs();
  UpdateAdjacentNNIsAfterAddAcceptedNNIs();
  // (6) Reset Accepted NNIs and GraftDAG.
  ClearAcceptedNNIs();
  RemoveAllGraftedNNIsFromDAG();
}

// ** NNI Evaluation

double NNIEngine::EvaluateNNI(const NNIOperation &nni) {
  double score = ComputeNNILikelihood(nni);
  return score;
}

double NNIEngine::ComputeNNILikelihood(const NNIOperation &nni) {
  Failwith("Currently no implementation for NNIEngine::ComputeNNILikelihood().");
  return 0.0f;
}

void NNIEngine::EvaluateAdjacentNNIs() {
  // Resize engine.
  ResizeGPEngineForGraftDAG();
  for (const auto &nni : GetAdjacentNNIs()) {
    const auto score = EvaluateNNI(nni);
    AddScoreForNNI(nni, score);
  }
}

// ** Add NNIs / Scores

void NNIEngine::GraftAdjacentNNIsToDAG() {
  for (const auto &nni : GetAdjacentNNIs()) {
    graft_dag_->AddGraftNodePair(nni);
  }
}

void NNIEngine::RemoveAllGraftedNNIsFromDAG() { graft_dag_->RemoveAllGrafts(); }

void NNIEngine::AddScoreForNNI(const NNIOperation &nni, const double score) {
  SafeInsert(scored_nnis_, nni, score);
};

double NNIEngine::GetScoreForNNI(const NNIOperation &nni) const {
  Assert(scored_nnis_.find(nni) != scored_nnis_.end(),
         "NNIEngine::GetScoreForNNI(): Score does not exist for requested NNI.");
  double score = scored_nnis_.at(nni);
  return score;
}

// ** Filters

void NNIEngine::FilterInit() {
  Failwith("Currently no implementation for NNIEngine::FilterInit().");
  // !ISSUE #405: Need to implement filtering scheme
}

void NNIEngine::FilterUpdate() {
  Failwith("Currently no implementation for NNIEngine::FilterUpdate().");
  // !ISSUE #405: Need to implement filtering scheme
}

void NNIEngine::FilterAdjacentNNIsToAcceptedNNIs() {
  Failwith(
      "Currently no implementation for NNIEngine::FilterAdjacentNNIsToAcceptedNNIs().");
  // !ISSUE #405: Need to implement filtering scheme
}

// ** DAG & GPEngine Maintenance

void NNIEngine::InitGPEngine() {
  Failwith("Currently no implementation for NNIEngine::InitGPEngine().");
  // !ISSUE #404: Need to be able to resize engine.
  gp_engine_.ProcessOperations(dag_.PopulatePLVs());
  gp_engine_.ProcessOperations(dag_.ComputeLikelihoods());
}

void NNIEngine::ResizeGPEngineForDAG() {
  Failwith("Currently no implementation for NNIEngine::ResizeGPEngineForDAG().");
  // !ISSUE #404: Need to be able to resize engine.
}

void NNIEngine::ResizeGPEngineForGraftDAG() {
  Failwith("Currently no implementation for NNIEngine::ResizeGPEngineForDAG().");
  // !ISSUE #404: Need to be able to resize engine.
}

void NNIEngine::PrepGPEngineForLikelihoods() {
  gp_engine_.ProcessOperations(dag_.PopulatePLVs());
  gp_engine_.ProcessOperations(dag_.ComputeLikelihoods());
}

void NNIEngine::AddAcceptedNNIsToDAG() {
  for (const auto &nni : accepted_nnis_) {
    dag_.AddNodePair(nni);
  }
}

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
  const size_t parent_id = dag_.GetDAGNodeId(parent_bitset);
  const size_t child_id = dag_.GetDAGNodeId(child_bitset);
  // Every new edge added is a potential new NNI.
  // Iterate over the parent and child node of the new pair.
  for (const size_t &node_id : {parent_id, child_id}) {
    // Get nodes adjacent to current node from both left and right edges.
    for (const bool is_edge_leafward : {true, false}) {
      // Get nodes adjacent to current node from both leafward and rootward directions.
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

void NNIEngine::UpdateAdjacentNNIsAfterAddAcceptedNNIs() {
  for (const auto &nni : GetAcceptedNNIs()) {
    const auto is_edge_on_left =
        Bitset::SubsplitIsChildOfWhichParentClade(nni.parent_, nni.child_);
    SafeAddOutputNNIsToAdjacentNNIs(nni.parent_, nni.child_, is_edge_on_left);
  }
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
    // If DAG already contains output parent and child nodes, and an edge between them,
    // then don't add it to the adjacent_nnis.
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

void NNIEngine::ClearAdjacentNNIs() {
  adjacent_nnis_.clear();
  scored_nnis_.clear();
}

void NNIEngine::ClearAcceptedNNIs() { accepted_nnis_.clear(); }
