// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "nni_evaluation_engine.hpp"

using KeyIndex = NNIEngine::KeyIndex;
using KeyIndexMap = NNIEngine::KeyIndexMap;
using KeyIndexMapPair = NNIEngine::KeyIndexMapPair;

// ** NNIEvaluationEngine

explicit NNIEvaluationEngine::NNIEvaluationEngine(NNIEngine &nni_engine)
    : nni_engine_(&nni_engine),
      dag_(&nni_engine.GetGPDAG()),
      graft_dag_(&nni_engine.GetGraftDAG()) {}

// ** NNIEvaluationEngineViaGP

NNIEvaluationEngineViaGP::NNIEvaluationEngineViaGP(NNIEngine &nni_engine,
                                                   GPEngine &gp_engine)
    : NNIEvaluationEngine(nni_engine), gp_engine_(&gp_engine) {}

void NNIEvaluationEngineViaGP::Init() {
  GetGPEngine().GrowPLVs(GetDAG().NodeCountWithoutDAGRoot());
  GetGPEngine().GrowGPCSPs(GetDAG().EdgeCountWithLeafSubsplits());
  GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());
  GetGPEngine().ProcessOperations(GetDAG().ComputeLikelihoods());
}

void NNIEvaluationEngineViaGP::Prep() {
  GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());
  GetGPEngine().ProcessOperations(GetDAG().ComputeLikelihoods());
}

void NNIEvaluationEngineViaGP::GrowEngineForDAG(
    std::optional<Reindexer> node_reindexer, std::optional<Reindexer> edge_reindexer) {
  GetGPEngine().GrowPLVs(GetDAG().NodeCountWithoutDAGRoot(), node_reindexer);
  GetGPEngine().GrowGPCSPs(GetDAG().EdgeCountWithLeafSubsplits(), edge_reindexer);
}

void NNIEvaluationEngineViaGP::GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                                         const bool via_reference,
                                                         const bool use_unique_temps) {
  if (via_reference) {
    if (use_unique_temps) {
      GetGPEngine().GrowSparePLVs(2 * adjacent_nnis.size());
    } else {
      GetGPEngine().GrowSparePLVs(2);
    }
    GetGPEngine().GrowSpareGPCSPs(adjacent_nnis.size());
  } else {
    GetGPEngine().GrowPLVs(GetGraftDAG().NodeCountWithoutDAGRoot());
    GetGPEngine().GrowGPCSPs(GetGraftDAG().EdgeCountWithLeafSubsplits());
  }
}

void NNIEvaluationEngineViaGP::ScoreAdjacentNNIs(const NNISet &adjacent_nnis) {
  // Compute Likelihoods
  const GPOperationVector operations =
      BuildGPOperationsForAdjacentNNILikelihoods(adjacent_nnis, true);
  GetGPEngine().ProcessOperations(operations);
  // Retrieve results from GPEngine and store in Scored NNIs.
  const auto proposed_nni_likelihoods =
      GetGPEngine().GetSparePerGPCSPLogLikelihoods(0, adjacent_nnis.size());
  size_t nni_count = 0;
  for (const auto &nni : adjacent_nnis) {
    GetScoredNNIs()[nni] = proposed_nni_likelihoods[nni_count];
    nni_count++;
  }
}

GPOperationVector NNIEvaluationEngineViaGP::BuildGPOperationsForAdjacentNNILikelihoods(
    const NNISet &adjacent_nnis, const bool via_reference) {
  GPOperationVector operations = GPOperationVector();
  GrowGPEngineForAdjacentNNILikelihoods(adjacent_nnis, via_reference);
  size_t nni_count = 0;
  for (const auto &nni : adjacent_nnis) {
    const auto pre_nni = GetDAG().FindNNINeighborInDAG(nni);
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

GPOperationVector NNIEvaluationEngineViaGP::BuildGPOperationsForNNILikelihood(
    const NNIOperation &nni, const NNIEngine::KeyIndexMap &nni_key_idx) const {
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

GPOperationVector NNIEvaluationEngineViaGP::BuildGPOperationsForAdjacentNNILikelihoods(
    const NNISet &adjacent_nnis, const bool via_reference) {
  GPOperationVector operations = GPOperationVector();
  GrowGPEngineForAdjacentNNILikelihoods(adjacent_nnis, via_reference);
  size_t nni_count = 0;
  for (const auto &nni : adjacent_nnis) {
    const auto pre_nni = GetDAG().FindNNINeighborInDAG(nni);
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

void NNIEvaluationEngineViaGP::GrowGPEngineForAdjacentNNILikelihoods(
    const NNISet &adjacent_nnis, const bool via_reference,
    const bool use_unique_temps) {
  if (via_reference) {
    if (use_unique_temps) {
      GetGPEngine().GrowSparePLVs(2 * adjacent_nnis.size());
    } else {
      GetGPEngine().GrowSparePLVs(2);
    }
    GetGPEngine().GrowSpareGPCSPs(adjacent_nnis.size());
  } else {
    GetGPEngine().GrowPLVs(GetGraftDAG().NodeCountWithoutDAGRoot());
    GetGPEngine().GrowGPCSPs(GetGraftDAG().EdgeCountWithLeafSubsplits());
  }
}

KeyIndexMapPair NNIEvaluationEngineViaGP::PassGPEngineDataFromPreNNIToPostNNIViaCopy(
    const NNIOperation &pre_nni, const NNIOperation &post_nni) {
  // Find data in pre-NNI.
  const auto pre_key_idx =
      GetNNIEngine().BuildKeyIndexMapForNNI(pre_nni, GetGraftDAG().NodeCount() - 1);
  // Find data in NNI.
  const auto &post_key_idx =
      GetNNIEngine().BuildKeyIndexMapForNNI(post_nni, GetGraftDAG().NodeCount() - 1);

  // Array for mapping from pre-NNI plvs to post-NNI plvs.
  const auto key_type_map =
      GetNNIEngine().BuildKeyIndexTypePairsFromPreNNIToPostNNI(pre_nni, post_nni);
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

KeyIndexMapPair
NNIEvaluationEngineViaGP::PassGPEngineDataFromPreNNIToPostNNIViaReference(
    const NNIOperation &pre_nni, const NNIOperation &post_nni, const size_t nni_count,
    const bool use_unique_temps) {
  // Find data in pre-NNI.
  KeyIndexMap pre_key_idx = GetNNIEngine().BuildKeyIndexMapForNNI(
      pre_nni, GetDAG().NodeCountWithoutDAGRoot());
  // Find data in post-NNI.
  KeyIndexMap post_key_idx =
      GetNNIEngine().BuildKeyIndexMapForPostNNIViaReferencePreNNI(pre_nni, post_nni,
                                                                  pre_key_idx);

  // Assign temporary place in engine to store intermediate values.
  size_t temp_offset_1 = 0;
  size_t temp_offset_2 = 0;
  if (use_unique_temps) {
    temp_offset_1 = nni_count * 2;
    temp_offset_2 = (nni_count * 2) + 1;
  } else {
    temp_offset_1 = 0;
    temp_offset_2 = 1;
  }
  post_key_idx[KeyIndex::Parent_RFocal] =
      GetGPEngine().GetPLVHandler().GetSparePVIndex(PVId(temp_offset_1)).value_;
  post_key_idx[KeyIndex::Child_P] =
      GetGPEngine().GetPLVHandler().GetSparePVIndex(PVId(temp_offset_2)).value_;
  post_key_idx[KeyIndex::Edge] = GetGPEngine().GetSpareGPCSPIndex(nni_count);

  // Copy over central edge data.
  for (const auto idx_type : {KeyIndex::Edge}) {
    const auto pre_edge_idx = EdgeId(pre_key_idx[idx_type]);
    const auto post_edge_idx = EdgeId(post_key_idx[idx_type]);
    GetGPEngine().CopyGPCSPData(pre_edge_idx, post_edge_idx);
  }

  return {pre_key_idx, post_key_idx};
}

// ** NNIEvaluationEngineViaTP

void NNIEvaluationEngineViaTP::Init() {
  GetTPEngine().GrowNodeData(GetDAG().NodeCount());
  GetTPEngine().GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits());
}

void NNIEvaluationEngineViaTP::Prep() {}

void NNIEvaluationEngineViaTP::GrowEngineForDAG(
    std::optional<Reindexer> node_reindexer, std::optional<Reindexer> edge_reindexer) {
  GetTPEngine().GrowNodeData(GetDAG().NodeCount(), node_reindexer);
  GetTPEngine().GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits(), edge_reindexer);
}

void NNIEvaluationEngineViaTP::GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                                         const bool via_reference,
                                                         const bool use_unique_temps) {
  if (via_reference) {
    if (use_unique_temps) {
      GetTPEngine().GrowSpareNodeData(2 * adjacent_nnis.size());
    } else {
      GetTPEngine().GrowSpareNodeData(2);
    }
    GetTPEngine().GrowSpareEdgeData(adjacent_nnis.size());
  } else {
    GetTPEngine().GrowNodeData(GetGraftDAG().NodeCountWithoutDAGRoot());
    GetTPEngine().GrowEdgeData(GetGraftDAG().EdgeCountWithLeafSubsplits());
  }
}

void NNIEvaluationEngineViaTP::ScoreAdjacentNNIs(const NNISet &adjacent_nnis) {
  // Retrieve results from TPEngine and store in Scored NNIs.
  for (const auto &nni : adjacent_nnis) {
    const auto pre_nni = GetDAG().FindNNINeighborInDAG(nni);
    GetScoredNNIs()[nni] =
        GetTPEngine().GetTopTreeLikelihoodWithProposedNNI(nni, pre_nni);
  }
}
