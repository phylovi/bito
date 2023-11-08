// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "nni_evaluation_engine.hpp"
#include "nni_engine.hpp"
#include "tp_engine.hpp"
#include "gp_engine.hpp"

// ** NNIEvalEngine

NNIEvalEngine::NNIEvalEngine(NNIEngine &nni_engine)
    : nni_engine_(&nni_engine),
      dag_(&nni_engine.GetDAG()),
      graft_dag_(&nni_engine.GetGraftDAG()) {}

double NNIEvalEngine::GetScoreByNNI(const NNIOperation &nni) const {
  auto it = GetScoredNNIs().find(nni);
  Assert(it != GetScoredNNIs().end(), "NNI does not exist in NNI Evaluation Engine.");
  return it->second;
}

double NNIEvalEngine::GetScoreByEdge(const EdgeId edge_id) const {
  auto nni = GetDAG().GetNNI(edge_id);
  return GetScoreByNNI(nni);
}

double NNIEvalEngine::GetMaxScore() const {
  double max = -INFINITY;
  for (const auto &[nni, score] : GetScoredNNIs()) {
    std::ignore = nni;
    if (max < score) {
      max = score;
    }
  }
  return max;
}

double NNIEvalEngine::GetMinScore() const {
  double min = INFINITY;
  for (const auto &[nni, score] : GetScoredNNIs()) {
    std::ignore = nni;
    if (min > score) {
      min = score;
    }
  }
  return min;
}

// ** NNIEvalEngineViaGP

NNIEvalEngineViaGP::NNIEvalEngineViaGP(NNIEngine &nni_engine, GPEngine &gp_engine)
    : NNIEvalEngine(nni_engine), gp_engine_(&gp_engine) {}

void NNIEvalEngineViaGP::Init() {
  GetGPEngine().GrowPLVs(GetDAG().NodeCountWithoutDAGRoot());
  GetGPEngine().GrowGPCSPs(GetDAG().EdgeCountWithLeafSubsplits());
}

void NNIEvalEngineViaGP::Prep() {
  GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());
  GetGPEngine().ProcessOperations(GetDAG().ComputeLikelihoods());
}

void NNIEvalEngineViaGP::GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                                          std::optional<Reindexer> edge_reindexer) {
  // Remove DAGRoot from node reindexing (for GPEngine).
  const Reindexer node_reindexer_without_root =
      node_reindexer.value().RemoveNewIndex(GetDAG().GetDAGRootNodeId().value_);
  GetGPEngine().GrowPLVs(GetDAG().NodeCountWithoutDAGRoot(),
                         node_reindexer_without_root);
  GetGPEngine().GrowGPCSPs(GetDAG().EdgeCountWithLeafSubsplits(), edge_reindexer);
}

void NNIEvalEngineViaGP::GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                                   const bool via_reference,
                                                   const bool use_unique_temps) {
  if (via_reference) {
    if (use_unique_temps) {
      GetGPEngine().GrowSparePLVs(GetSpareNodesPerNNI() * adjacent_nnis.size());
    } else {
      GetGPEngine().GrowSparePLVs(GetSpareNodesPerNNI());
    }
    GetGPEngine().GrowSpareGPCSPs(adjacent_nnis.size());
  } else {
    GetGPEngine().GrowPLVs(GetGraftDAG().NodeCountWithoutDAGRoot());
    GetGPEngine().GrowGPCSPs(GetGraftDAG().EdgeCountWithLeafSubsplits());
  }
}

void NNIEvalEngineViaGP::UpdateEngineAfterModifyingDAG(
    const std::map<NNIOperation, NNIOperation> &pre_nni_to_nni,
    const size_t prev_node_count, const Reindexer &node_reindexer,
    const size_t prev_edge_count, const Reindexer &edge_reindexer) {
  using namespace GPOperations;
  const bool copy_branch_lengths = false;

  auto &branch_handler = GetGPEngine().GetBranchLengthHandler();
  // Find all new edge ids.
  std::set<EdgeId> new_edge_ids;
  for (size_t i = prev_edge_count; i < edge_reindexer.size(); i++) {
    EdgeId edge_id = EdgeId(edge_reindexer.GetNewIndexByOldIndex(i));
    new_edge_ids.insert(edge_id);
  }

  // Grow and Reindex GPEngine.
  GetGPEngine().GrowPLVs(node_reindexer.size());
  GetGPEngine().GrowGPCSPs(edge_reindexer.size());
  for (size_t i = prev_edge_count; i < edge_reindexer.size(); i++) {
    const EdgeId edge_id = edge_reindexer.GetOldIndexByNewIndex(i);
    branch_handler(edge_id) = branch_handler.GetDefaultBranchLength();
  }
  // Copy over branch lengths from pre-NNI to post-NNI.
  if (copy_branch_lengths) {
    for (const auto &[pre_nni, nni] : pre_nni_to_nni) {
      CopyGPEngineDataAfterAddingNNI(pre_nni, nni);
    }
  }
  // Update SBN Priors.
  auto sbn_prior = GetDAG().BuildUniformOnTopologicalSupportPrior();
  auto unconditional_node_probabilities =
      GetDAG().UnconditionalNodeProbabilities(sbn_prior);
  auto inverted_sbn_prior =
      GetDAG().InvertedGPCSPProbabilities(sbn_prior, unconditional_node_probabilities);
  GetGPEngine().InitializePriors(std::move(sbn_prior),
                                 std::move(unconditional_node_probabilities.segment(
                                     0, GetDAG().NodeCountWithoutDAGRoot())),
                                 std::move(inverted_sbn_prior));
  if (use_null_priors_) {
    GetGPEngine().SetNullPrior();
  }
  // Update PLVs.
  GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());

  // Optimize branch lengths.
  if (IsOptimizeNewEdges()) {
    for (const auto &[pre_nni, nni] : pre_nni_to_nni) {
      std::ignore = pre_nni;
      NNIBranchLengthOptimization(nni, new_edge_ids);
    }
    GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());
  }

  GetGPEngine().ProcessOperations(GetDAG().ComputeLikelihoods());
  auto likelihoods = GetGPEngine().GetPerGPCSPLogLikelihoods();
}

void NNIEvalEngineViaGP::CopyGPEngineDataAfterAddingNNI(const NNIOperation &pre_nni,
                                                        const NNIOperation &nni) {
  auto &branch_handler = GetGPEngine().GetBranchLengthHandler();
  auto CopyBranchLengths = [this, &branch_handler](
                               const auto pre_node_id, const auto post_node_id,
                               const Direction dir, const SubsplitClade clade) {
    for (const auto adj_node_id :
         GetDAG().GetDAGNode(pre_node_id).GetNeighbors(dir, clade)) {
      const auto pre_parent_id =
          (dir == Direction::Rootward) ? adj_node_id : pre_node_id;
      const auto post_parent_id =
          (dir == Direction::Rootward) ? adj_node_id : post_node_id;
      const auto pre_child_id =
          (dir == Direction::Rootward) ? pre_node_id : adj_node_id;
      const auto post_child_id =
          (dir == Direction::Rootward) ? post_node_id : adj_node_id;
      const auto pre_edge_id = GetDAG().GetEdgeIdx(pre_parent_id, pre_child_id);
      const auto post_edge_id = GetDAG().GetEdgeIdx(post_parent_id, post_child_id);
      GetGPEngine().CopyGPCSPData(pre_edge_id, post_edge_id);
    }
  };

  // Build mapping according to the NNI swap.
  const auto post_parent_id = GetDAG().GetDAGNodeId(nni.GetParent());
  const auto post_child_id = GetDAG().GetDAGNodeId(nni.GetChild());
  const auto post_edge_id = GetDAG().GetEdgeIdx(post_parent_id, post_child_id);
  const auto pre_parent_id = GetDAG().GetDAGNodeId(pre_nni.GetParent());
  const auto pre_child_id = GetDAG().GetDAGNodeId(pre_nni.GetChild());
  const auto pre_edge_id = GetDAG().GetEdgeIdx(pre_parent_id, pre_child_id);
  const auto pre_sister_clade = GetDAG().GetSisterClade(pre_edge_id);
  const auto clade_map = NNIOperation::BuildNNICladeMapFromPreNNIToNNI(nni, pre_nni);

  NNICladeEnum::Array<NodeId> post_node_ids;
  post_node_ids[clade_map[NNIClade::ParentFocal]] = post_parent_id;
  post_node_ids[clade_map[NNIClade::ParentSister]] = post_parent_id;
  post_node_ids[clade_map[NNIClade::ChildLeft]] = post_child_id;
  post_node_ids[clade_map[NNIClade::ChildRight]] = post_child_id;
  DoubleVector parents, sisters, children;
  // Copy parent edges.
  for (const auto clade : SubsplitCladeEnum::Iterator()) {
    CopyBranchLengths(pre_parent_id, post_node_ids[NNIClade::ParentFocal],
                      Direction::Rootward, clade);
  }
  // Copy sister edges.
  CopyBranchLengths(pre_parent_id, post_node_ids[NNIClade::ParentSister],
                    Direction::Leafward, pre_sister_clade);
  // Copy central edge.
  branch_handler.Get(post_edge_id) = branch_handler.Get(pre_edge_id);
  // Copy left and right child edges.
  for (const auto clade : SubsplitCladeEnum::Iterator()) {
    const auto nni_clade =
        (clade == SubsplitClade::Left) ? NNIClade::ChildLeft : NNIClade::ChildRight;
    CopyBranchLengths(pre_child_id, post_node_ids[nni_clade], Direction::Leafward,
                      clade);
  }
}

void NNIEvalEngineViaGP::ScoreAdjacentNNIs(const NNISet &adjacent_nnis) {
  ComputeAdjacentNNILikelihoods(adjacent_nnis, true);
}

double NNIEvalEngineViaGP::ScoreInternalNNIByNNI(const NNIOperation &nni) const {
  Assert(GetDAG().ContainsNNI(nni), "DAG does not contain NNI.");
  const auto edge_id = GetDAG().GetEdgeIdx(nni);
  return ScoreInternalNNIByEdge(edge_id);
}

double NNIEvalEngineViaGP::ScoreInternalNNIByEdge(const EdgeId &edge_id) const {
  return GetGPEngine().GetPerGPCSPLogLikelihoods(edge_id.value_, 1)[0];
}

size_t NNIEvalEngineViaGP::GetSpareNodesPerNNI() const { return spare_nodes_per_nni_; }

size_t NNIEvalEngineViaGP::GetSpareEdgesPerNNI() const { return spare_edges_per_nni_; }

void NNIEvalEngineViaGP::ComputeAdjacentNNILikelihoods(const NNISet &adjacent_nnis,
                                                       const bool via_reference) {
  GrowGPEngineForAdjacentNNILikelihoods(adjacent_nnis, via_reference);
  std::vector<double> scored_nnis(adjacent_nnis.size());
  size_t nni_count = 0;
  size_t offset = 0;
  for (const auto &nni : adjacent_nnis) {
    auto [likelihood, new_offset] = ComputeAdjacentNNILikelihood(nni, offset);
    std::ignore = likelihood;
    std::ignore = new_offset;
    nni_count++;
  }
}

std::pair<double, size_t> NNIEvalEngineViaGP::ComputeAdjacentNNILikelihood(
    const NNIOperation &nni, const size_t offset) {
  using namespace GPOperations;
  GPOperationVector ops;
  auto &pvs = GetGPEngine().GetPLVHandler();
  const NNIOperation pre_nni = GetDAG().FindNNINeighborInDAG(nni);

  // const std::set<NodeId> rootsplit_node_ids{GetDAG().GetRootsplitNodeIds().begin(),
  //                                           GetDAG().GetRootsplitNodeIds().end()};
  std::set<NodeId> rootsplit_node_ids;
  for (const auto node_id : GetDAG().GetRootsplitNodeIds()) {
    rootsplit_node_ids.insert(node_id);
  }

  // Get temp PVs for adjacent NNI.
  const auto pv_ids = GetTempAdjPVIds();
  // Get mapped nodes and temp edges for adjacent NNI.
  auto [node_ids, edge_ids] =
      GetMappedAdjNodeIdsAndTempAdjEdgeIds(pre_nni, nni, IsCopyNewEdges());
  size_t edge_count = edge_ids.parents.size() + edge_ids.sisters.size() +
                      edge_ids.leftchildren.size() + edge_ids.rightchildren.size() + 1;
  size_t new_offset = edge_count;
  GetGPEngine().GrowSpareGPCSPs(edge_count);

  // Rootward Pass
  auto UpdateLeftChildRootward = [&]() {
    // Evolve up leftchild edge: leftchild_p -> child_phatleft.
    size_t i = 0;
    ops.push_back(ZeroPLV{pv_ids.child_phatleft.value_});
    for (const auto adj_node_id : node_ids.grand_leftchildren) {
      const PVId leftchild_p = pvs.GetPVIndex(PLVType::P, adj_node_id);
      ops.push_back(IncrementWithWeightedEvolvedPLV{pv_ids.child_phatleft.value_,
                                                    edge_ids.leftchildren[i].value_,
                                                    leftchild_p.value_});
      i++;
    }
  };
  auto UpdateRightChildRootward = [&]() {
    // Evolve up rightchild edge: rightchild_p -> child_phatright.
    size_t i = 0;
    ops.push_back(ZeroPLV{pv_ids.child_phatright.value_});
    for (const auto adj_node_id : node_ids.grand_rightchildren) {
      const PVId rightchild_p = pvs.GetPVIndex(PLVType::P, adj_node_id);
      ops.push_back(IncrementWithWeightedEvolvedPLV{pv_ids.child_phatright.value_,
                                                    edge_ids.rightchildren[i].value_,
                                                    rightchild_p.value_});
      i++;
    }
  };
  auto UpdateCentralRootward = [&]() {
    // child_p = child_phatleft \circ child_phatright..
    ops.push_back(Multiply{pv_ids.child_p.value_, pv_ids.child_phatleft.value_,
                           pv_ids.child_phatright.value_});
    // Evolve up central edge: child_p -> parent_phatfocal.
    ops.push_back(ZeroPLV{pv_ids.parent_phatfocal.value_});
    ops.push_back(IncrementWithWeightedEvolvedPLV{pv_ids.parent_phatfocal.value_,
                                                  edge_ids.central.value_,
                                                  pv_ids.child_p.value_});
  };
  auto UpdateSisterRootward = [&]() {
    // Evolve up sister edge: sister_p -> parent_phatsister.
    size_t i = 0;
    ops.push_back(ZeroPLV{pv_ids.parent_phatsister.value_});
    for (const auto adj_node_id : node_ids.grand_sisters) {
      const PVId sister_p = pvs.GetPVIndex(PLVType::P, adj_node_id);
      ops.push_back(IncrementWithWeightedEvolvedPLV{pv_ids.parent_phatsister.value_,
                                                    edge_ids.sisters[i].value_,
                                                    sister_p.value_});
      i++;
    }
  };
  auto UpdateParentRootward = [&]() {
    // parent_p = parent_phatfocal \circ parent_phatsister.
    ops.push_back(Multiply{pv_ids.parent_p.value_, pv_ids.parent_phatfocal.value_,
                           pv_ids.parent_phatsister.value_});
  };
  auto NNIRootwardPass = [&]() {
    UpdateLeftChildRootward();
    UpdateRightChildRootward();
    UpdateCentralRootward();
    UpdateSisterRootward();
    UpdateParentRootward();
  };
  // Leafward Pass
  auto UpdateLeftChildLeafward = [&]() {
    // child_rleft = child_rhat \circ child_rhatright.
    ops.push_back(Multiply{pv_ids.child_rleft.value_, pv_ids.child_rhat.value_,
                           pv_ids.child_phatright.value_});
  };
  auto UpdateRightChildLeafward = [&]() {
    // child_rright = child_rhat \circ child_rhatleft.
    ops.push_back(Multiply{pv_ids.child_rright.value_, pv_ids.child_rhat.value_,
                           pv_ids.child_phatleft.value_});
  };
  auto UpdateCentralLeafward = [&]() {
    // parent_rfocal = parent_rhat \circ parent_rhatsister.
    ops.push_back(Multiply{pv_ids.parent_rfocal.value_, pv_ids.parent_rhat.value_,
                           pv_ids.parent_phatsister.value_});
    // Evolve down central edge: parent_rfocal -> child_rhat.
    ops.push_back(ZeroPLV{pv_ids.child_rhat.value_});
    ops.push_back(IncrementWithWeightedEvolvedPLV{pv_ids.child_rhat.value_,
                                                  edge_ids.central.value_,
                                                  pv_ids.parent_rfocal.value_});
  };
  auto UpdateSisterLeafward = [&]() {
    // parent_rsister = parent_rhat \circ parent_rhatfocal.
    ops.push_back(Multiply{pv_ids.parent_rsister.value_, pv_ids.parent_rhat.value_,
                           pv_ids.parent_phatfocal.value_});
  };
  auto UpdateParentLeafward = [&]() {
    // Evolve down parent edge: grandparent_rfocal -> parent_rhat.
    ops.push_back(ZeroPLV{pv_ids.parent_rhat.value_});
    size_t i = 0;
    const bool is_dag_root = GetDAG().IsNodeRoot(node_ids.grand_parents[0]);
    if (is_dag_root) {
      ops.push_back(SetToStationaryDistribution{pv_ids.parent_rhat.value_,
                                                edge_ids.parents[0].value_});
    } else {
      for (const auto adj_node_id : node_ids.grand_parents) {
        const EdgeId edge_id = GetDAG().GetEdgeIdx(adj_node_id, node_ids.parent_focal);
        const auto focal_clade = GetDAG().GetFocalClade(edge_id);
        const PVId grandparent_rfocal =
            pvs.GetPVIndex(PLVTypeEnum::RPLVType(focal_clade), adj_node_id);
        ops.push_back(IncrementWithWeightedEvolvedPLV{pv_ids.parent_rhat.value_,
                                                      edge_ids.parents[i].value_,
                                                      grandparent_rfocal.value_});
        i++;
      }
    }
  };
  auto NNILeafwardPass = [&]() {
    UpdateParentLeafward();
    UpdateCentralLeafward();
    UpdateSisterLeafward();
    UpdateLeftChildLeafward();
    UpdateRightChildLeafward();
  };
  // Branch Length Optimization
  auto OptimizeLeftChild = [&](bool do_update = true) {
    size_t i = 0;
    for (const auto adj_node_id : node_ids.grand_leftchildren) {
      const PVId grand_leftchild_p = pvs.GetPVIndex(PLVType::P, adj_node_id);
      ops.push_back(OptimizeBranchLength{pv_ids.child_rleft.value_,
                                         grand_leftchild_p.value_,
                                         edge_ids.leftchildren[i].value_});
      i++;
    }
    if (do_update) UpdateLeftChildRootward();
  };
  auto OptimizeRightChild = [&](bool do_update = true) {
    size_t i = 0;
    for (const auto adj_node_id : node_ids.grand_rightchildren) {
      const PVId grand_rightchild_p = pvs.GetPVIndex(PLVType::P, adj_node_id);
      ops.push_back(OptimizeBranchLength{pv_ids.child_rright.value_,
                                         grand_rightchild_p.value_,
                                         edge_ids.rightchildren[i].value_});
      i++;
    }
    if (do_update) UpdateRightChildRootward();
  };
  auto OptimizeCentral = [&](bool do_update = true) {
    if (do_update) UpdateCentralLeafward();
    ops.push_back(OptimizeBranchLength{pv_ids.parent_rfocal.value_,
                                       pv_ids.child_p.value_, edge_ids.central.value_});
    if (do_update) UpdateCentralRootward();
  };
  auto OptimizeSister = [&](bool do_update = true) {
    if (do_update) UpdateSisterLeafward();
    size_t i = 0;
    for (const auto adj_node_id : node_ids.grand_sisters) {
      const PVId grand_sister_p = pvs.GetPVIndex(PLVType::P, adj_node_id);
      ops.push_back(OptimizeBranchLength{pv_ids.parent_rsister.value_,
                                         grand_sister_p.value_,
                                         edge_ids.sisters[i].value_});
      i++;
    }
    if (do_update) UpdateSisterRootward();
  };
  auto OptimizeParent = [&](bool do_update = true) {
    if (do_update) UpdateParentLeafward();
    size_t i = 0;
    if ((node_ids.grand_parents.size() == 1) &&
        GetDAG().IsNodeRoot(node_ids.grand_parents[0])) {
      return;
    } else {
      for (const auto adj_node_id : node_ids.grand_parents) {
        const EdgeId edge_id = GetDAG().GetEdgeIdx(adj_node_id, node_ids.parent_focal);
        const auto focal_clade = GetDAG().GetFocalClade(edge_id);
        const PVId grandparent_rfocal =
            pvs.GetPVIndex(PLVTypeEnum::RPLVType(focal_clade), adj_node_id);
        ops.push_back(OptimizeBranchLength{
            grandparent_rfocal.value_,
            pv_ids.parent_p.value_,
            edge_ids.parents[i].value_,
        });
        i++;
      }
    }
    if (do_update) UpdateParentRootward();
  };
  auto NNIBranchLengthOptimization = [&](bool do_update = true) {
    OptimizeLeftChild(do_update);
    OptimizeRightChild(do_update);
    OptimizeSister(do_update);
    OptimizeCentral(do_update);
    OptimizeParent(do_update);
  };

  // Branch Length Optimization.
  if (IsOptimizeNewEdges()) {
    NNIRootwardPass();
    NNILeafwardPass();
    GetGPEngine().ProcessOperations(ops);
    ops.clear();

    NNIBranchLengthOptimization();
    NNILeafwardPass();
    for (size_t iter = 0; iter < GetOptimizationMaxIteration(); iter++) {
      GetGPEngine().ProcessOperations(ops);
    }
    ops.clear();
  }

  // Compute likelihood of central edge.
  NNIRootwardPass();
  NNILeafwardPass();
  ops.push_back(GPOperations::Likelihood{
      edge_ids.central.value_, pv_ids.parent_rfocal.value_, pv_ids.child_p.value_});
  GetGPEngine().ProcessOperations(ops);
  ops.clear();

  double likelihood =
      GetGPEngine().GetPerGPCSPLogLikelihoods(edge_ids.central.value_, 1)[0];
  GetScoredNNIs()[nni] = likelihood;
  return {likelihood, new_offset};
}

NNIEvalEngineViaGP::AdjNodeAndEdgeIds NNIEvalEngineViaGP::GetAdjNodeAndEdgeIds(
    const NNIOperation &nni) const {
  AdjNodeIds node_ids;
  AdjEdgeIds edge_ids;
  Assert(GetDAG().ContainsNNI(nni), "Given NNI does not exist in DAG.");

  const auto AddAdjNodesAndEdgesToVector =
      [this](NodeIdVector &adj_node_ids, EdgeIdVector &adj_edge_ids,
             const NodeId node_id, const Direction dir, const SubsplitClade clade) {
        for (const auto adj_node_id :
             GetDAG().GetDAGNode(node_id).GetNeighbors(dir, clade)) {
          if (adj_node_id >= GetDAG().NodeCount()) return;
          const auto parent_id = (dir == Direction::Rootward) ? adj_node_id : node_id;
          const auto child_id = (dir == Direction::Rootward) ? node_id : adj_node_id;
          const auto edge_id = GetDAG().GetEdgeIdx(parent_id, child_id);
          adj_node_ids.push_back(adj_node_id);
          adj_edge_ids.push_back(edge_id);
        }
      };

  node_ids.parent_focal = GetDAG().GetDAGNodeId(nni.GetParent());
  node_ids.parent_sister = node_ids.parent_focal;
  node_ids.child_left = GetDAG().GetDAGNodeId(nni.GetChild());
  node_ids.child_right = node_ids.child_left;
  edge_ids.central = GetDAG().GetEdgeIdx(node_ids.parent_focal, node_ids.child_left);
  const auto sister_clade = GetDAG().GetSisterClade(edge_ids.central);

  for (const auto clade : SubsplitCladeEnum::Iterator()) {
    AddAdjNodesAndEdgesToVector(node_ids.grand_parents, edge_ids.parents,
                                node_ids.parent_focal, Direction::Rootward, clade);
  }
  AddAdjNodesAndEdgesToVector(node_ids.grand_sisters, edge_ids.sisters,
                              node_ids.parent_sister, Direction::Leafward,
                              sister_clade);
  AddAdjNodesAndEdgesToVector(node_ids.grand_leftchildren, edge_ids.leftchildren,
                              node_ids.child_left, Direction::Leafward,
                              SubsplitClade::Left);
  AddAdjNodesAndEdgesToVector(node_ids.grand_rightchildren, edge_ids.rightchildren,
                              node_ids.child_right, Direction::Leafward,
                              SubsplitClade::Right);
  return {node_ids, edge_ids};
}

NNIEvalEngineViaGP::AdjNodeAndEdgeIds
NNIEvalEngineViaGP::GetMappedAdjNodeIdsAndTempAdjEdgeIds(
    const NNIOperation &pre_nni, const NNIOperation &nni,
    const bool copy_branch_lengths) {
  AdjNodeIds node_ids;
  // Get edges and nodes from pre-NNI.
  auto [pre_node_ids, pre_edge_ids] = GetAdjNodeAndEdgeIds(pre_nni);
  const auto clade_map = NNIOperation::BuildNNICladeMapFromPreNNIToNNI(nni, pre_nni);

  // Remap adj node ids from pre-NNI to post-NNI.
  NNICladeEnum::Array<NodeIdVector *> unmapped_node_ids, mapped_node_ids;
  unmapped_node_ids[NNIClade::ParentFocal] = &pre_node_ids.grand_parents;
  unmapped_node_ids[NNIClade::ParentSister] = &pre_node_ids.grand_sisters;
  unmapped_node_ids[NNIClade::ChildLeft] = &pre_node_ids.grand_leftchildren;
  unmapped_node_ids[NNIClade::ChildRight] = &pre_node_ids.grand_rightchildren;
  for (const auto nni_clade : NNICladeEnum::Iterator()) {
    mapped_node_ids[nni_clade] = unmapped_node_ids[clade_map[nni_clade]];
  }
  node_ids.grand_parents = *mapped_node_ids[NNIClade::ParentFocal];
  node_ids.grand_sisters = *mapped_node_ids[NNIClade::ParentSister];
  node_ids.grand_leftchildren = *mapped_node_ids[NNIClade::ChildLeft];
  node_ids.grand_rightchildren = *mapped_node_ids[NNIClade::ChildRight];

  // Remap central node ids from pre-NNI to post-NNI.
  NNICladeEnum::Array<NodeId> unmapped_node_id, mapped_node_id;
  unmapped_node_id[NNIClade::ParentFocal] = pre_node_ids.parent_focal;
  unmapped_node_id[NNIClade::ParentSister] = pre_node_ids.parent_sister;
  unmapped_node_id[NNIClade::ChildLeft] = pre_node_ids.child_left;
  unmapped_node_id[NNIClade::ChildRight] = pre_node_ids.child_right;
  for (const auto nni_clade : NNICladeEnum::Iterator()) {
    mapped_node_id[nni_clade] = unmapped_node_id[clade_map[nni_clade]];
  }
  node_ids.parent_focal = mapped_node_id[NNIClade::ParentFocal];
  node_ids.parent_sister = mapped_node_id[NNIClade::ParentSister];
  node_ids.child_left = mapped_node_id[NNIClade::ChildLeft];
  node_ids.child_right = mapped_node_id[NNIClade::ChildRight];

  // Remap edges from pre-NNI to post-NNI.
  AdjEdgeIds edge_ids;
  NNICladeEnum::Array<EdgeIdVector *> unmapped_edge_ids, mapped_edge_ids;
  unmapped_edge_ids[NNIClade::ParentFocal] = &pre_edge_ids.parents;
  unmapped_edge_ids[NNIClade::ParentSister] = &pre_edge_ids.sisters;
  unmapped_edge_ids[NNIClade::ChildLeft] = &pre_edge_ids.leftchildren;
  unmapped_edge_ids[NNIClade::ChildRight] = &pre_edge_ids.rightchildren;
  for (const auto nni_clade : NNICladeEnum::Iterator()) {
    mapped_edge_ids[nni_clade] = unmapped_edge_ids[clade_map[nni_clade]];
  }
  edge_ids.central = pre_edge_ids.central;
  edge_ids.parents = *mapped_edge_ids[NNIClade::ParentFocal];
  edge_ids.sisters = *mapped_edge_ids[NNIClade::ParentSister];
  edge_ids.leftchildren = *mapped_edge_ids[NNIClade::ChildLeft];
  edge_ids.rightchildren = *mapped_edge_ids[NNIClade::ChildRight];

  // Assign temp edges for post-NNI.
  AdjEdgeIds temp_edge_ids = GetTempAdjEdgeIds(node_ids);

  // Copy pre-NNI branches to post-NNI branches.
  auto CopyEdgeData = [&](const EdgeId pre_edge_id, const EdgeId post_edge_id) {
    GetGPEngine().CopyGPCSPData(pre_edge_id, post_edge_id);
  };
  auto CopyEdgesData = [&](const EdgeIdVector pre_edge_ids,
                           const EdgeIdVector post_edge_ids) {
    for (size_t i = 0; i < pre_edge_ids.size(); i++) {
      CopyEdgeData(pre_edge_ids[i], post_edge_ids[i]);
    }
  };

  if (copy_branch_lengths) {
    CopyEdgeData(edge_ids.central, temp_edge_ids.central);
    CopyEdgesData(edge_ids.parents, temp_edge_ids.parents);
    CopyEdgesData(edge_ids.sisters, temp_edge_ids.sisters);
    CopyEdgesData(edge_ids.leftchildren, temp_edge_ids.leftchildren);
    CopyEdgesData(edge_ids.rightchildren, temp_edge_ids.rightchildren);
  }

  return {node_ids, temp_edge_ids};
}

NNIEvalEngineViaGP::AdjEdgeIds NNIEvalEngineViaGP::GetTempAdjEdgeIds(
    const AdjNodeIds &node_ids) {
  AdjEdgeIds edge_ids;
  // Grow spare edges if neccessary.
  const size_t spare_edge_count =
      1 + node_ids.grand_parents.size() + node_ids.grand_sisters.size() +
      node_ids.grand_leftchildren.size() + node_ids.grand_rightchildren.size();
  GetGPEngine().GrowSpareGPCSPs(spare_edge_count);

  // Assign temp edges based on number of adjacent edges in pre-NNI.
  size_t edge_offset = 0;
  const auto GetNextSpareEdgeIds = [this, &edge_offset](
                                       const NodeIdVector &adj_node_ids,
                                       const NodeId node_id, const Direction dir) {
    EdgeIdVector edge_ids;
    for (size_t i = 0; i < adj_node_ids.size(); i++) {
      const EdgeId edge_id = GetGPEngine().GetSpareGPCSPIndex(edge_offset);
      edge_offset++;
      edge_ids.push_back(edge_id);
    }
    return edge_ids;
  };

  // Assign central edge to first temp index.
  edge_ids.central = GetGPEngine().GetSpareGPCSPIndex(edge_offset);
  edge_offset++;
  // Assign adjacent edges.
  edge_ids.parents = GetNextSpareEdgeIds(node_ids.grand_parents, node_ids.parent_focal,
                                         Direction::Rootward);
  edge_ids.sisters = GetNextSpareEdgeIds(node_ids.grand_sisters, node_ids.parent_sister,
                                         Direction::Leafward);
  edge_ids.leftchildren = GetNextSpareEdgeIds(node_ids.grand_leftchildren,
                                              node_ids.child_left, Direction::Leafward);
  edge_ids.rightchildren = GetNextSpareEdgeIds(
      node_ids.grand_rightchildren, node_ids.child_right, Direction::Leafward);

  return edge_ids;
}

NNIEvalEngineViaGP::AdjPVIds NNIEvalEngineViaGP::GetTempAdjPVIds() {
  AdjPVIds pv_ids;
  // Grow PLV handler if necessary.
  GetGPEngine().GrowSparePLVs(2 * PLVTypeEnum::Count);
  // Get next available spare PVId.
  size_t node_offset = 0;
  const auto GetNextSparePVId = [this, &node_offset]() {
    auto &pvs = GetGPEngine().GetPLVHandler();
    const PVId pvid = pvs.GetSparePVIndex(PVId(node_offset));
    node_offset++;
    return pvid;
  };
  // Get PLVs for parent.
  pv_ids.parent_p = GetNextSparePVId();
  pv_ids.parent_phatfocal = GetNextSparePVId();
  pv_ids.parent_phatsister = GetNextSparePVId();
  pv_ids.parent_rhat = GetNextSparePVId();
  pv_ids.parent_rfocal = GetNextSparePVId();
  pv_ids.parent_rsister = GetNextSparePVId();
  // Get PLVs for child.
  pv_ids.child_p = GetNextSparePVId();
  pv_ids.child_phatleft = GetNextSparePVId();
  pv_ids.child_phatright = GetNextSparePVId();
  pv_ids.child_rhat = GetNextSparePVId();
  pv_ids.child_rleft = GetNextSparePVId();
  pv_ids.child_rright = GetNextSparePVId();
  return pv_ids;
}

std::set<EdgeId> NNIEvalEngineViaGP::BuildSetOfEdgeIdsAdjacentToNNI(
    const NNIOperation &nni) const {
  std::set<EdgeId> edge_ids;
  const auto parent_id = GetDAG().GetDAGNodeId(nni.GetParent());
  const auto child_id = GetDAG().GetDAGNodeId(nni.GetChild());
  const auto central_edge_id = GetDAG().GetEdgeIdx(parent_id, child_id);
  const auto sister_clade = GetDAG().GetSisterClade(central_edge_id);
  // Collect left children.
  for (const auto adj_node_id : GetDAG().GetDAGNode(child_id).GetNeighbors(
           Direction::Leafward, SubsplitClade::Left)) {
    const auto edge_id = GetDAG().GetEdgeIdx(child_id, adj_node_id);
    edge_ids.insert(edge_id);
  }
  // Collect right children.
  for (const auto adj_node_id : GetDAG().GetDAGNode(child_id).GetNeighbors(
           Direction::Leafward, SubsplitClade::Right)) {
    const auto edge_id = GetDAG().GetEdgeIdx(child_id, adj_node_id);
    edge_ids.insert(edge_id);
  }
  // Collect sisters.
  for (const auto adj_node_id :
       GetDAG().GetDAGNode(parent_id).GetNeighbors(Direction::Leafward, sister_clade)) {
    const auto edge_id = GetDAG().GetEdgeIdx(parent_id, adj_node_id);
    edge_ids.insert(edge_id);
  }
  // Collect central.
  edge_ids.insert(central_edge_id);
  // Collect parents.
  for (const auto clade : SubsplitCladeEnum::Iterator()) {
    for (const auto adj_node_id :
         GetDAG().GetDAGNode(parent_id).GetNeighbors(Direction::Rootward, clade)) {
      if (GetDAG().IsNodeRoot(adj_node_id)) {
        continue;
      }
      const auto edge_id = GetDAG().GetEdgeIdx(adj_node_id, parent_id);
      edge_ids.insert(edge_id);
    }
  }
  return edge_ids;
}

std::set<Bitset> NNIEvalEngineViaGP::BuildSetOfPCSPsAdjacentToNNI(
    const NNIOperation &nni) const {
  std::set<Bitset> pcsps;
  auto edge_ids = BuildSetOfEdgeIdsAdjacentToNNI(nni);
  for (const auto edge_id : edge_ids) {
    pcsps.insert(GetDAG().GetDAGEdgeBitset(edge_id));
  }
  return pcsps;
}

void NNIEvalEngineViaGP::GrowGPEngineForAdjacentNNILikelihoods(
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

NNIEvalEngine::KeyIndexMapPair
NNIEvalEngineViaGP::PassGPEngineDataFromPreNNIToPostNNIViaCopy(
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

NNIEvalEngine::KeyIndexMapPair
NNIEvalEngineViaGP::PassGPEngineDataFromPreNNIToPostNNIViaReference(
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

void NNIEvalEngineViaGP::BranchLengthOptimization() {
  const auto ops = GetDAG().BranchLengthOptimization();
  for (size_t iter = 0; iter < GetOptimizationMaxIteration(); iter++) {
    GetGPEngine().ProcessOperations(ops);
  }
}

void NNIEvalEngineViaGP::BranchLengthOptimization(
    const std::set<EdgeId> &edges_to_optimize) {
  const auto ops = GetDAG().BranchLengthOptimization(edges_to_optimize);
  for (size_t iter = 0; iter < GetOptimizationMaxIteration(); iter++) {
    GetGPEngine().ProcessOperations(ops);
    GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());
  }
}

void NNIEvalEngineViaGP::NNIBranchLengthOptimization(
    const NNIOperation &nni, const std::set<EdgeId> &new_edge_ids) {
  using namespace GPOperations;
  const auto adj_edge_ids = BuildSetOfEdgeIdsAdjacentToNNI(nni);
  // Get the edges that are new to the DAG.
  std::set<EdgeId> new_adj_edge_ids;
  for (const auto edge_id : adj_edge_ids) {
    if (new_edge_ids.find(edge_id) != new_edge_ids.end()) {
      new_adj_edge_ids.insert(edge_id);
    }
  }
  // Initial optimization.
  GPOperationVector init_ops;
  auto BuildOpForRootward = [this, &init_ops](const NodeId node_id) {
    auto node = GetDAG().GetDAGNode(node_id);
    auto &pvs = GetGPEngine().GetPLVHandler();
    if (GetDAG().IsNodeLeaf(node_id)) return;
    // Update PHatLeft.
    init_ops.push_back(ZeroPLV{pvs.GetPVIndex(PLVType::PHatLeft, node_id).value_});
    for (const auto adj_node_id :
         node.GetNeighbors(Direction::Leafward, SubsplitClade::Left)) {
      const auto edge = GetDAG().GetDAGEdge(GetDAG().GetEdgeIdx(node_id, adj_node_id));
      init_ops.push_back(IncrementWithWeightedEvolvedPLV{
          pvs.GetPVIndex(PLVType::PHatLeft, node_id).value_, edge.GetId().value_,
          pvs.GetPVIndex(PLVType::P, adj_node_id).value_});
    }
    // Update PHatRight.
    init_ops.push_back(ZeroPLV{pvs.GetPVIndex(PLVType::PHatRight, node_id).value_});
    for (const auto adj_node_id :
         node.GetNeighbors(Direction::Leafward, SubsplitClade::Right)) {
      const auto edge = GetDAG().GetDAGEdge(GetDAG().GetEdgeIdx(node_id, adj_node_id));
      init_ops.push_back(IncrementWithWeightedEvolvedPLV{
          pvs.GetPVIndex(PLVType::PHatRight, node_id).value_, edge.GetId().value_,
          pvs.GetPVIndex(PLVType::P, adj_node_id).value_});
    }
    // Update P.
    init_ops.push_back(Multiply{pvs.GetPVIndex(PLVType::P, node_id).value_,
                                pvs.GetPVIndex(PLVType::PHatRight, node_id).value_,
                                pvs.GetPVIndex(PLVType::PHatLeft, node_id).value_});
  };
  auto BuildOpForLeafward = [this, &init_ops](const NodeId node_id) {
    auto node = GetDAG().GetDAGNode(node_id);
    auto &pvs = GetGPEngine().GetPLVHandler();
    if (GetDAG().IsNodeRoot(node_id)) return;
    // Update RHat.
    init_ops.push_back(ZeroPLV{pvs.GetPVIndex(PLVType::RHat, node_id).value_});
    for (const auto clade : SubsplitCladeEnum::Iterator()) {
      for (const auto adj_node_id : node.GetNeighbors(Direction::Rootward, clade)) {
        const auto edge =
            GetDAG().GetDAGEdge(GetDAG().GetEdgeIdx(adj_node_id, node_id));
        init_ops.push_back(IncrementWithWeightedEvolvedPLV{
            pvs.GetPVIndex(PLVType::RHat, node_id).value_, edge.GetId().value_,
            pvs.GetPVIndex(PLVTypeEnum::RPLVType(edge.GetSubsplitClade()), adj_node_id)
                .value_});
      }
    }
    // Update RLeft.
    init_ops.push_back(Multiply{pvs.GetPVIndex(PLVType::RRight, node_id).value_,
                                pvs.GetPVIndex(PLVType::RHat, node_id).value_,
                                pvs.GetPVIndex(PLVType::PHatLeft, node_id).value_});
    // Update RRight.
    init_ops.push_back(Multiply{pvs.GetPVIndex(PLVType::RLeft, node_id).value_,
                                pvs.GetPVIndex(PLVType::RHat, node_id).value_,
                                pvs.GetPVIndex(PLVType::PHatRight, node_id).value_});
  };
  auto BuildOpForBranchLengthOptimization =
      [this, &init_ops, &new_adj_edge_ids, &BuildOpForRootward, &BuildOpForLeafward](
          const EdgeId edge_id, const bool do_update = true) {
        if (GetDAG().IsEdgeRoot(edge_id)) return;
        if (new_adj_edge_ids.find(edge_id) == new_adj_edge_ids.end()) return;
        const auto &edge = GetDAG().GetDAGEdge(edge_id);
        auto &pvs = GetGPEngine().GetPLVHandler();
        if (do_update) {
          BuildOpForRootward(edge.GetParent());
          BuildOpForLeafward(edge.GetParent());
        }
        PVId parent_rfocal = pvs.GetPVIndex(
            PLVTypeEnum::RPLVType(edge.GetSubsplitClade()), edge.GetParent());
        PVId child_p = pvs.GetPVIndex(PLVType::P, edge.GetChild());
        init_ops.push_back(
            OptimizeBranchLength{parent_rfocal.value_, child_p.value_, edge_id.value_});
        if (do_update) {
          BuildOpForRootward(edge.GetParent());
        }
      };
  const auto [node_ids, edge_ids] = GetAdjNodeAndEdgeIds(nni);
  std::ignore = node_ids;
  // Build Optimization Operations.
  for (const auto edge_id : edge_ids.leftchildren) {
    BuildOpForBranchLengthOptimization(edge_id);
  }
  for (const auto edge_id : edge_ids.rightchildren) {
    BuildOpForBranchLengthOptimization(edge_id);
  }
  for (const auto edge_id : edge_ids.sisters) {
    BuildOpForBranchLengthOptimization(edge_id);
  }
  BuildOpForBranchLengthOptimization(edge_ids.central);
  for (const auto edge_id : edge_ids.parents) {
    BuildOpForBranchLengthOptimization(edge_id);
  }
  // Build update PLVs.
  // Rootward Pass.
  for (const auto edge_id : edge_ids.leftchildren) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForRootward(edge.GetChild());
  }
  for (const auto edge_id : edge_ids.rightchildren) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForRootward(edge.GetChild());
  }
  for (const auto edge_id : edge_ids.sisters) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForRootward(edge.GetChild());
  }
  for (const auto edge_id : {edge_ids.central}) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForRootward(edge.GetChild());
  }
  for (const auto edge_id : edge_ids.parents) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForRootward(edge.GetChild());
  }
  // Leafward Pass.
  for (const auto edge_id : edge_ids.parents) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForLeafward(edge.GetParent());
  }
  for (const auto edge_id : {edge_ids.central}) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForLeafward(edge.GetParent());
  }
  for (const auto edge_id : edge_ids.sisters) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForLeafward(edge.GetParent());
  }
  for (const auto edge_id : edge_ids.rightchildren) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForLeafward(edge.GetParent());
  }
  for (const auto edge_id : edge_ids.leftchildren) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    BuildOpForLeafward(edge.GetParent());
  }

  for (size_t iter = 0; iter < GetOptimizationMaxIteration(); iter++) {
    GetGPEngine().ProcessOperations(init_ops);
    GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());
  }
}

void NNIEvalEngineViaGP::NNIBranchLengthOptimization(const NNIOperation &nni) {
  using namespace GPOperations;
  const auto adj_edge_ids = BuildSetOfEdgeIdsAdjacentToNNI(nni);
  // Initial optimization.
  GPOperationVector init_ops;
  auto BuildGPOpForBranchLengthOptimization = [this, &init_ops](const EdgeId edge_id) {
    if (GetDAG().IsEdgeRoot(edge_id)) return;
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    auto &pvs = GetGPEngine().GetPLVHandler();
    PVId parent_rfocal = pvs.GetPVIndex(PLVTypeEnum::RPLVType(edge.GetSubsplitClade()),
                                        edge.GetParent());
    PVId child_p = pvs.GetPVIndex(PLVType::P, edge.GetChild());
    init_ops.push_back(
        OptimizeBranchLength{parent_rfocal.value_, child_p.value_, edge_id.value_});
  };
  const auto [node_ids, edge_ids] = GetAdjNodeAndEdgeIds(nni);
  std::ignore = node_ids;
  for (const auto edge_id : edge_ids.leftchildren) {
    BuildGPOpForBranchLengthOptimization(edge_id);
  }
  for (const auto edge_id : edge_ids.rightchildren) {
    BuildGPOpForBranchLengthOptimization(edge_id);
  }
  for (const auto edge_id : edge_ids.sisters) {
    BuildGPOpForBranchLengthOptimization(edge_id);
  }
  BuildGPOpForBranchLengthOptimization(edge_ids.central);
  for (const auto edge_id : edge_ids.parents) {
    BuildGPOpForBranchLengthOptimization(edge_id);
  }
  GetGPEngine().ProcessOperations(init_ops);
  // Optimize after updating entire DAG.
  const auto full_ops = GetDAG().BranchLengthOptimization(adj_edge_ids);
  for (size_t iter = 1; iter < GetOptimizationMaxIteration(); iter++) {
    GetGPEngine().ProcessOperations(full_ops);
    GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());
  }
}

// ** NNIEvalEngineViaTP

void NNIEvalEngineViaTP::Init() { GetTPEngine().Initialize(); }

void NNIEvalEngineViaTP::Prep() {
  GetTPEngine().InitializeChoiceMap();
  GetTPEngine().InitializeScores();
}

void NNIEvalEngineViaTP::GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                                          std::optional<Reindexer> edge_reindexer) {
  GetTPEngine().GrowNodeData(GetDAG().NodeCount(), node_reindexer);
  GetTPEngine().GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits(), edge_reindexer);
}

void NNIEvalEngineViaTP::GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                                   const bool via_reference,
                                                   const bool use_unique_temps) {
  if (via_reference) {
    if (use_unique_temps) {
      GetTPEngine().GrowSpareNodeData(GetSpareNodesPerNNI() * adjacent_nnis.size());
      GetTPEngine().GrowSpareEdgeData(GetSpareEdgesPerNNI() * adjacent_nnis.size());
    } else {
      GetTPEngine().GrowSpareNodeData(GetSpareNodesPerNNI());
      GetTPEngine().GrowSpareEdgeData(GetSpareEdgesPerNNI());
    }
  } else {
    GetTPEngine().GrowNodeData(GetGraftDAG().NodeCountWithoutDAGRoot());
    GetTPEngine().GrowEdgeData(GetGraftDAG().EdgeCountWithLeafSubsplits());
  }
}

void NNIEvalEngineViaTP::UpdateEngineAfterModifyingDAG(
    const std::map<NNIOperation, NNIOperation> &pre_nni_to_nni,
    const size_t prev_node_count, const Reindexer &node_reindexer,
    const size_t prev_edge_count, const Reindexer &edge_reindexer) {
  GetTPEngine().UpdateAfterModifyingDAG(pre_nni_to_nni, prev_node_count, node_reindexer,
                                        prev_edge_count, edge_reindexer);
}

void NNIEvalEngineViaTP::ScoreAdjacentNNIs(const NNISet &adjacent_nnis) {
  const auto best_edge_map =
      GetTPEngine().BuildMapOfProposedNNIPCSPsToBestPreNNIEdges(adjacent_nnis);
  for (const auto &nni : adjacent_nnis) {
    const auto pre_nni = GetDAG().FindNNINeighborInDAG(nni);
    GetScoredNNIs()[nni] =
        GetTPEngine().GetTopTreeScoreWithProposedNNI(nni, pre_nni, 0, best_edge_map);
  }
}

double NNIEvalEngineViaTP::ScoreInternalNNIByNNI(const NNIOperation &nni) const {
  Assert(GetDAG().ContainsNNI(nni), "DAG does not contain NNI.");
  const auto edge_id = GetDAG().GetEdgeIdx(nni);
  return ScoreInternalNNIByEdge(edge_id);
}

double NNIEvalEngineViaTP::ScoreInternalNNIByEdge(const EdgeId &edge_id) const {
  return GetTPEngine().GetTopTreeScore(edge_id);
}

size_t NNIEvalEngineViaTP::GetSpareNodesPerNNI() const {
  return GetTPEngine().GetSpareNodesPerNNI();
}

size_t NNIEvalEngineViaTP::GetSpareEdgesPerNNI() const {
  return GetTPEngine().GetSpareNodesPerNNI();
}
