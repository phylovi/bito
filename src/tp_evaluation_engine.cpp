// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "tp_evaluation_engine.hpp"
#include "tp_engine.hpp"
#include "pv_handler.hpp"
#include "numerical_utils.hpp"

// ** TPEvalEngine

TPEvalEngine::TPEvalEngine(TPEngine &tp_engine)
    : tp_engine_(&tp_engine),
      dag_(&GetTPEngine().GetDAG()),
      graft_dag_(&GetTPEngine().GetGraftDAG()),
      site_pattern_(&GetTPEngine().GetSitePattern()) {
  GrowNodeData(GetDAG().NodeCount(), std::nullopt, std::nullopt, true);
  GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits(), std::nullopt, std::nullopt, true);
}

void TPEvalEngine::Initialize() { Failwith("Pure virtual function call."); }

void TPEvalEngine::ComputeScores(std::optional<EdgeIdVector> opt_edge_ids) {
  Failwith("Pure virtual function call.");
}

void TPEvalEngine::GrowNodeData(const size_t new_node_count,
                                std::optional<const Reindexer> node_reindexer,
                                std::optional<const size_t> explicit_alloc,
                                const bool on_init) {
  // No node data to grow/reindex.
}

void TPEvalEngine::GrowEdgeData(const size_t edge_count,
                                std::optional<const Reindexer> edge_reindexer,
                                std::optional<const size_t> explicit_alloc,
                                const bool on_init) {
  // Build resizer for resizing data.
  Resizer resizer(GetTPEngine().GetEdgeCount(), GetTPEngine().GetSpareEdgeCount(),
                  GetTPEngine().GetAllocatedEdgeCount(), edge_count, std::nullopt,
                  explicit_alloc, GetTPEngine().GetResizingFactor());
  resizer.ApplyResizeToEigenVector<EigenVectorXd, double>(GetTopTreeScores(),
                                                          DOUBLE_NEG_INF);
  // Reindex work space to realign with DAG.
  if (edge_reindexer.has_value()) {
    Reindexer::ReindexInPlace<EigenVectorXd, double>(
        top_tree_per_edge_, edge_reindexer.value(), resizer.GetNewCount());
  }
}

void TPEvalEngine::GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                                    std::optional<Reindexer> edge_reindexer) {
  GrowNodeData(GetDAG().NodeCount(), node_reindexer, std::nullopt, false);
  GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits(), edge_reindexer, std::nullopt,
               false);
}

void TPEvalEngine::GrowSpareNodeData(const size_t new_node_spare_count) {
  if (new_node_spare_count > GetTPEngine().GetSpareNodeCount()) {
    GetTPEngine().GrowSpareNodeData(new_node_spare_count);
  }
  GrowNodeData(GetTPEngine().GetNodeCount());
}

void TPEvalEngine::GrowSpareEdgeData(const size_t new_edge_spare_count) {
  if (new_edge_spare_count > GetTPEngine().GetSpareEdgeCount()) {
    GetTPEngine().GrowSpareEdgeData(new_edge_spare_count);
  }
  GrowEdgeData(GetTPEngine().GetEdgeCount());
}

void TPEvalEngine::GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                             const bool via_reference,
                                             const bool use_unique_temps) {
  Failwith("Pure virtual function call.");
}

void TPEvalEngine::UpdateEngineAfterDAGAddNodePair(const NNIOperation &post_nni,
                                                   const NNIOperation &pre_nni,
                                                   std::optional<size_t> new_tree_id) {
  Failwith("Pure virtual function call.");
}

void TPEvalEngine::UpdateEngineAfterModifyingDAG(
    const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
    const size_t prev_node_count, const Reindexer &node_reindexer,
    const size_t prev_edge_count, const Reindexer &edge_reindexer) {
  Failwith("Pure virtual function.");
}

double TPEvalEngine::GetTopTreeScoreWithEdge(const EdgeId edge_id) const {
  double score = top_tree_per_edge_[edge_id.value_];
  return score;
}

double TPEvalEngine::GetTopTreeScoreWithProposedNNI(
    const NNIOperation &post_nni, const NNIOperation &pre_nni,
    const size_t spare_offset, std::optional<BitsetEdgeIdMap> best_edge_map) {
  Failwith("Pure virtual function call.");
  return DOUBLE_NEG_INF;
}

void TPEvalEngine::CopyEdgeData(const EdgeId src_edge_id, const EdgeId dest_edge_id) {
  top_tree_per_edge_[dest_edge_id.value_] = top_tree_per_edge_[src_edge_id.value_];
}

// ** TPEvalEngineViaLikelihood

TPEvalEngineViaLikelihood::TPEvalEngineViaLikelihood(TPEngine &tp_engine,
                                                     const std::string &mmap_path)
    : TPEvalEngine(tp_engine),
      likelihood_pvs_(mmap_path, GetDAG().EdgeCountWithLeafSubsplits(),
                      GetSitePattern().PatternCount(), 2.0),
      branch_handler_(tp_engine.GetDAG()) {
  GrowNodeData(GetDAG().NodeCount(), std::nullopt, std::nullopt, true);
  GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits(), std::nullopt, std::nullopt, true);
  InitializeBranchLengthHandler();
}

void TPEvalEngineViaLikelihood::Initialize() {
  // Set all PVs to Zero
  ZeroPVs();
  // Populate Leaves with Site Patterns.
  PopulateLeafPVsWithSitePatterns();
  // Populate Rootsplit with Stationary Distribution.
  PopulateRootPVsWithStationaryDistribution();
  // Populate rootward and leafward PVs.
  PopulatePVs();
}

void TPEvalEngineViaLikelihood::ZeroPVs() {
  for (EdgeId edge_id = 0; edge_id < GetPVs().GetCount(); edge_id++) {
    for (const auto pv_type : PLVTypeEnum::Iterator()) {
      GetPVs().GetPV(pv_type, edge_id).setZero();
    }
  }
}

void TPEvalEngineViaLikelihood::PopulatePVs() {
  // Rootward Pass (populate P-PVs)
  PopulateRootwardPVs();
  // Leafward Pass (populate R-PVs)
  PopulateLeafwardPVs();
}

void TPEvalEngineViaLikelihood::PopulateRootwardPVs() {
  const auto rootward_node_ids = GetDAG().RootwardNodeTraversalTrace(false);
  for (const auto node_id : rootward_node_ids) {
    PopulateRootwardPVForNode(node_id);
  }

  // TODO instead of using nodes, switch to edges?
  // const auto rootward_edge_ids = GetDAG().RootwardEdgeTraversalTrace(true);
  // for (const auto edge_id : rootward_edge_ids) {
  //   PopulateRootwardPVForEdge(edge_id);
  // }
}

void TPEvalEngineViaLikelihood::PopulateLeafwardPVs() {
  const auto leafward_node_ids = GetDAG().LeafwardNodeTraversalTrace(true);
  for (const auto node_id : leafward_node_ids) {
    PopulateLeafwardPVForNode(node_id);
  }

  // TODO instead of using nodes, switch to edges?
  // const auto leafward_edge_ids = GetDAG().LeafwardEdgeTraversalTrace(true);
  // for (const auto edge_id : leafward_edge_ids) {
  //   PopulateLeafwardPVForEdge(edge_id);
  // }
}

void TPEvalEngineViaLikelihood::GrowNodeData(
    const size_t node_count, std::optional<const Reindexer> node_reindexer,
    std::optional<const size_t> explicit_alloc, const bool on_init) {
  // No node data to resize.
}

void TPEvalEngineViaLikelihood::GrowEdgeData(
    const size_t edge_count, std::optional<const Reindexer> edge_reindexer,
    std::optional<const size_t> explicit_alloc, const bool on_init) {
  bool is_quiet = true;
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cout);
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);
  // Build resizer for resizing data.
  Resizer resizer(GetTPEngine().GetEdgeCount(), GetTPEngine().GetSpareEdgeCount(),
                  GetTPEngine().GetAllocatedEdgeCount(), edge_count, std::nullopt,
                  explicit_alloc, GetTPEngine().GetResizingFactor());
  GetDAGBranchHandler().Resize(resizer.GetNewCount(), edge_reindexer);
  GetMatrix().conservativeResize(resizer.GetNewAlloc(),
                                 GetTPEngine().GetSitePattern().PatternCount());
  GetMatrix().conservativeResize(resizer.GetNewPadded(),
                                 GetTPEngine().GetSitePattern().PatternCount());
  resizer.ApplyResizeToEigenVector<EigenVectorXd, double>(GetTopTreeScores(),
                                                          DOUBLE_NEG_INF);
  GetPVs().Resize(resizer.GetNewCount(), resizer.GetNewAlloc(), resizer.GetNewSpare());

  // Reindex work space to realign with DAG.
  if (edge_reindexer.has_value()) {
    auto pv_reindexer = GetPVs().BuildPVReindexer(
        edge_reindexer.value(), resizer.GetOldCount(), resizer.GetNewCount());
    GetPVs().Reindex(pv_reindexer);
    os << "TPLikelihood::ReindexEdgeData::PVs: " << timer.Lap() << std::endl;
    Reindexer::ReindexInPlace<EigenVectorXd, double>(
        GetTopTreeScores(), edge_reindexer.value(), resizer.GetNewCount());
    os << "TPLikelihood::ReindexEdgeData::TopTreeScores: " << timer.Lap() << std::endl;
  }
}

void TPEvalEngineViaLikelihood::GrowSpareNodeData(const size_t new_node_spare_count) {
  if (new_node_spare_count > GetTPEngine().GetSpareNodeCount()) {
    GetTPEngine().GrowSpareNodeData(new_node_spare_count);
  }
  GrowNodeData(GetTPEngine().GetNodeCount());
}

void TPEvalEngineViaLikelihood::GrowSpareEdgeData(const size_t new_edge_spare_count) {
  if (new_edge_spare_count > GetTPEngine().GetSpareEdgeCount()) {
    GetTPEngine().GrowSpareEdgeData(new_edge_spare_count);
  }
  GrowEdgeData(GetTPEngine().GetEdgeCount());
}

void TPEvalEngineViaLikelihood::GrowEngineForDAG(
    std::optional<Reindexer> node_reindexer, std::optional<Reindexer> edge_reindexer) {
  const auto &dag = GetTPEngine().GetDAG();
  GrowNodeData(dag.NodeCount(), node_reindexer, std::nullopt, false);
  GrowEdgeData(dag.EdgeCountWithLeafSubsplits(), edge_reindexer, std::nullopt, false);
}

void TPEvalEngineViaLikelihood::GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                                          const bool via_reference,
                                                          const bool use_unique_temps) {
  if (via_reference) {
    if (use_unique_temps) {
      GrowSpareNodeData(spare_nodes_per_nni_ * adjacent_nnis.size());
      GrowSpareEdgeData(spare_edges_per_nni_ * adjacent_nnis.size());
    } else {
      GrowSpareNodeData(spare_nodes_per_nni_);
      GrowSpareEdgeData(spare_edges_per_nni_);
    }
    GrowSpareEdgeData(adjacent_nnis.size());
  } else {
    GrowNodeData(GetGraftDAG().NodeCountWithoutDAGRoot());
    GrowEdgeData(GetGraftDAG().EdgeCountWithLeafSubsplits());
  }
}

void TPEvalEngineViaLikelihood::UpdateEngineAfterDAGAddNodePair(
    const NNIOperation &post_nni, const NNIOperation &pre_nni,
    std::optional<size_t> new_tree_id) {
  // Copy over branch lengths.
  GetTPEngine().CopyOverEdgeDataFromPreNNIToPostNNI(
      post_nni, pre_nni,
      [this](const EdgeId src, const EdgeId dest) { CopyEdgeData(src, dest); },
      new_tree_id);

  const EdgeId focal_edge_id = GetDAG().GetEdgeIdx(post_nni);
  const auto &choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(focal_edge_id);
  // Center and adjacent edges, arranged rootward.
  const EdgeIdVector adj_edge_ids = {{choices.left_child, choices.right_child,
                                      focal_edge_id, choices.sister, choices.parent}};
  // Update PVs.
  for (const EdgeId edge_id : adj_edge_ids) {
    PopulateRootwardPVForEdge(edge_id);
  }
  for (auto it = adj_edge_ids.rbegin(); it != adj_edge_ids.rend(); ++it) {
    const EdgeId edge_id = *it;
    PopulateLeafwardPVForEdge(edge_id);
  }
  // Optimize focal branch and all adjacent branches.
  if (do_optimize_new_edges_) {
    for (const EdgeId edge_id : adj_edge_ids) {
      BranchLengthOptimization(edge_id, false);
    }
  }
}

bool watch_all = false;
std::set<std::string> nni_watchlist{{"0xFE0C0"}};
auto BitsetToHash = [](const Bitset &bitset) { return bitset.ToHashString(5); };

void TPEvalEngineViaLikelihood::UpdateEngineAfterModifyingDAG(
    const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
    const size_t prev_node_count, const Reindexer &node_reindexer,
    const size_t prev_edge_count, const Reindexer &edge_reindexer) {
  bool is_quiet = true;
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cout);
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);

  // Find edges to optimize. New edges were added to the DAG during current iteration.
  // NNI edges are are the central edge in NNIs. Extra edges are the new edges that are
  // not the choice edges of the NNIs. Update edges are edges that need their scores
  // recomputed.
  std::set<EdgeId> new_edges, nni_edges, extra_edges, update_edges;
  for (size_t i = prev_edge_count; i < edge_reindexer.size(); i++) {
    const EdgeId edge_id = EdgeId(edge_reindexer.GetNewIndexByOldIndex(i));
    const auto &choice = GetTPEngine().GetChoiceMap(edge_id);
    new_edges.insert(edge_id);
    extra_edges.insert(edge_id);
    update_edges.insert(edge_id);
    update_edges.insert(choice.parent);
    update_edges.insert(choice.sister);
    update_edges.insert(choice.left_child);
    update_edges.insert(choice.right_child);
  }
  for (const auto &[post_nni, pre_nni] : nni_to_pre_nni) {
    std::ignore = pre_nni;
    const auto edge_id = GetDAG().GetEdgeIdx(post_nni);
    const auto &choice = GetTPEngine().GetChoiceMap(edge_id);
    nni_edges.insert(edge_id);
    extra_edges.erase(choice.parent);
    extra_edges.erase(choice.sister);
    extra_edges.erase(edge_id);
    extra_edges.erase(choice.left_child);
    extra_edges.erase(choice.right_child);
    update_edges.insert(edge_id);
    update_edges.insert(choice.parent);
    update_edges.insert(choice.sister);
    update_edges.insert(choice.left_child);
    update_edges.insert(choice.right_child);
    update_edges.erase(EdgeId(NoId));
  }

  // Find topological sort of edges.
  EdgeIdVector rootward_edges(update_edges.begin(), update_edges.end());
  std::sort(rootward_edges.begin(), rootward_edges.end(),
            [this](const EdgeId lhs, const EdgeId rhs) {
              return GetDAG().GetDAGEdge(lhs).GetParent() <
                     GetDAG().GetDAGEdge(rhs).GetParent();
            });
  EdgeIdVector leafward_edges(update_edges.begin(), update_edges.end());
  std::sort(leafward_edges.begin(), leafward_edges.end(),
            [this](const EdgeId lhs, const EdgeId rhs) {
              return GetDAG().GetDAGEdge(lhs).GetChild() >
                     GetDAG().GetDAGEdge(rhs).GetChild();
            });

  // Populate Leaves with Site Patterns.
  PopulateLeafPVsWithSitePatterns();
  // Populate Rootsplit with Stationary Distribution.
  PopulateRootPVsWithStationaryDistribution();

  std::map<EdgeId, ProposedNNIInfo> nni_infos;
  for (const auto &[post_nni, pre_nni] : nni_to_pre_nni) {
    nni_infos[GetDAG().GetEdgeIdx(post_nni)] = GetRealNNIInfo(post_nni, pre_nni);
  }
  for (auto &[edge_id, nni_info] : nni_infos) {
    for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
      auto adj_edge_id = nni_info.adj_edge_ids[adj_nni_type];
      nni_info.do_optimize_edge[adj_nni_type] =
          (new_edges.find(adj_edge_id) != new_edges.end());
      // if (adj_edge_id != NoId) {
      //   if (nni_info.do_optimize_edge[adj_nni_type]) {
      //     auto ref_edge_id = nni_info.ref_edge_ids[adj_nni_type];
      //     branch_handler_(adj_edge_id) = branch_handler_(ref_edge_id);
      //   }
      // }
    }
  }

  os << "UpdateEngineAfterModifyingDAG::Init: " << timer.Lap() << std::endl;

  // bool is_on_watchlist = watch_all;
  // for (auto &[edge_id, nni_info] : nni_infos) {
  //   std::ignore = edge_id;
  //   is_on_watchlist |= nni_watchlist.find(nni_info.adj_pcsps.focal.ToHashString(5))
  //   !=
  //                      nni_watchlist.end();
  // }

  // // TODO Print Branch Lengths.
  // std::cout << "  === INIT_MODIFYING_DAG ===" << std::endl;
  // if (is_on_watchlist) {
  //   for (auto &[edge_id, nni_info] : nni_infos) {
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": PCSPS [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << nni_info.adj_pcsps[adj_nni_type].ToHashString(5) << " ";
  //     }

  //     std::cout << "]" << std::endl;
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": ADJ_CHOICES " <<
  //     std::endl; for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       auto edge_id = nni_info.adj_edge_ids[adj_nni_type];
  //       if (edge_id == NoId) {
  //         std::cout << "  [ ]" << std::endl;
  //         continue;
  //       }
  //       auto pcsps = GetTPEngine().GetChoiceMap().GetEdgeChoicePCSPs(edge_id);
  //       std::cout << "  [ ";
  //       for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //         std::cout << pcsps[adj_nni_type].ToHashString(5) << " ";
  //       }
  //       std::cout << "]" << std::endl;
  //     }

  //     std::cout << nni_info.post_nni.ToHashString(5) << ": REF_CHOICES " <<
  //     std::endl; for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << NNIAdjacentEnum::ToString(adj_nni_type) << " ";
  //       auto edge_id = nni_info.ref_edge_ids[adj_nni_type];
  //       if (edge_id == NoId) {
  //         std::cout << "  [ ]" << std::endl;
  //         continue;
  //       }
  //       auto pcsps = GetTPEngine().GetChoiceMap().GetEdgeChoicePCSPs(edge_id);
  //       std::cout << "  [ ";
  //       for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //         std::cout << pcsps[adj_nni_type].ToHashString(5) << " ";
  //       }
  //       std::cout << "]" << std::endl;
  //     }

  //     std::cout << nni_info.post_nni.ToHashString(5) << ": OLD_EDGES [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << !nni_info.do_optimize_edge[adj_nni_type] << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": PVS [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << GetPVs().ToHashString(nni_info.ref_primary_pv_ids[adj_nni_type],
  //       5)
  //                 << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": BEFORE_BLS [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       if (nni_info.adj_edge_ids[adj_nni_type] == NoId) continue;
  //       std::cout << DblToScientific(
  //                        branch_handler_(nni_info.adj_edge_ids[adj_nni_type]))
  //                 << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //   }
  // }

  // Initialize new PLVs.
  auto RootwardPass = [&]() {
    for (const auto edge_id : rootward_edges) {
      PopulateRootwardPVForEdge(edge_id);
    }
  };
  auto LeafwardPass = [&]() {
    for (const auto edge_id : leafward_edges) {
      PopulateLeafwardPVForEdge(edge_id);
    }
  };

  // // TODO can we remove this?
  // auto NNIRootwardPass = [&](const EdgeId edge_id) {
  //   const auto &choice = GetTPEngine().GetChoiceMap(edge_id);
  //   const auto pv_ids = GetLocalPVIdsOfEdge(edge_id);
  //   // Evolve up child P-PLVs.
  //   SetToEvolvedPV(pv_ids.child_phatleft_, choice.left_child, pv_ids.leftchild_p_);
  //   SetToEvolvedPV(pv_ids.child_phatright_, choice.right_child,
  //   pv_ids.rightchild_p_); MultiplyPVs(pv_ids.child_p_, pv_ids.child_phatleft_,
  //   pv_ids.child_phatright_);
  //   // Evolve up parent P-PLVs
  //   SetToEvolvedPV(pv_ids.parent_phatsister_, choice.sister, pv_ids.sister_p_);
  //   SetToEvolvedPV(pv_ids.parent_phatfocal_, edge_id, pv_ids.child_p_);
  //   MultiplyPVs(pv_ids.parent_p_, pv_ids.parent_phatfocal_,
  //   pv_ids.parent_phatsister_);
  // };
  // auto NNILeafwardPass = [&](const EdgeId edge_id) {
  //   const auto &choice = GetTPEngine().GetChoiceMap(edge_id);
  //   const auto pv_ids = GetLocalPVIdsOfEdge(edge_id);
  //   // If the parent is not the DAG root, then evolve grandparent down to parent.
  //   if (pv_ids.grandparent_rfocal_ != NoId) {
  //     SetToEvolvedPV(pv_ids.parent_rhat_, choice.parent, pv_ids.grandparent_rfocal_);
  //   }
  //   // Evolve down parent R-PLVs.
  //   MultiplyPVs(pv_ids.parent_rfocal_, pv_ids.parent_rhat_,
  //   pv_ids.parent_phatsister_); MultiplyPVs(pv_ids.parent_rsister_,
  //   pv_ids.parent_rhat_, pv_ids.parent_phatfocal_);
  //   SetToEvolvedPV(pv_ids.child_rhat_, edge_id, pv_ids.parent_rfocal_);
  //   // Update results.
  //   MultiplyPVs(pv_ids.child_rleft_, pv_ids.child_rhat_, pv_ids.child_phatright_);
  //   MultiplyPVs(pv_ids.child_rright_, pv_ids.child_rhat_, pv_ids.child_phatleft_);
  // };
  // auto NNIUpdatePVs = [&]() {
  //   // Update new NNI PVs.
  //   for (const auto edge_id : nni_edges) {
  //     NNIRootwardPass(edge_id);
  //     NNILeafwardPass(edge_id);
  //   }
  // };
  // std::ignore = NNIUpdatePVs;

  auto OptimizeEdge = [&](const NNIAdjacent nni_adj_type, const EdgeId edge_id,
                          const EdgeId parent_edge_id, const EdgeId focal_edge_id,
                          const bool is_not_child_edge = true,
                          const bool is_not_parent_edge = true,
                          const bool do_optimize = true) {
    const auto focal = GetDAG().GetFocalClade(edge_id);
    const auto sister = GetDAG().GetSisterClade(edge_id);
    if (is_not_child_edge) {
      // If child does not go to a leaf, update child_p (in case leftchild or
      // rightchild branch length were changed).
      MultiplyPVs(GetPVs().GetPVIndex(PLVType::P, edge_id),
                  GetPVs().GetPVIndex(PLVType::PHatLeft, edge_id),
                  GetPVs().GetPVIndex(PLVType::PHatRight, edge_id));
    }
    if (is_not_parent_edge) {
      // Update parent_rfocal(in case that sister branch length was changed).
      if (!GetDAG().IsEdgeRoot(edge_id)) {
        MultiplyPVs(GetPVs().GetPVIndex(PLVTypeEnum::RPLVType(focal), parent_edge_id),
                    GetPVs().GetPVIndex(PLVType::RHat, parent_edge_id),
                    GetPVs().GetPVIndex(PLVTypeEnum::PPLVType(sister), parent_edge_id));
      } else {
        TakePVValue(GetPVs().GetPVIndex(PLVTypeEnum::RPLVType(focal), parent_edge_id),
                    GetPVs().GetPVIndex(PLVType::RHat, parent_edge_id));
      }
    }
    // Optimize branch length.
    if (do_optimize) {
      const auto &[parent_rfocal, child_p] = GetPrimaryPVIdsOfEdge(edge_id);
      // if (is_on_watchlist) {
      //   std::cout << "DAG::optimize_edge::"
      //             << GetDAG().GetDAGEdgeBitset(edge_id).ToHashString(5) << " "
      //             << GetPVs().ToHashString(parent_rfocal, 5) << " "
      //             << GetPVs().ToHashString(child_p, 5) << std::endl;
      // }
      branch_handler_.OptimizeBranchLength(edge_id, parent_rfocal, child_p, false);
    }
    if (is_not_parent_edge) {
      // Update parent_phatfocal after changing branch length.
      SetToEvolvedPV(GetPVs().GetPVIndex(PLVTypeEnum::PPLVType(focal), parent_edge_id),
                     edge_id, GetPVs().GetPVIndex(PLVType::P, edge_id));
      // Update parent_p after changing branch length.
      MultiplyPVs(GetPVs().GetPVIndex(PLVType::P, parent_edge_id),
                  GetPVs().GetPVIndex(PLVType::PHatLeft, parent_edge_id),
                  GetPVs().GetPVIndex(PLVType::PHatRight, parent_edge_id));
    }
  };

  RootwardPass();
  LeafwardPass();

  // // TODO Print Branch Lengths.
  // std::cout << "  === BEFORE_MODIFYING_DAG ===" << std::endl;
  // if (is_on_watchlist) {
  //   for (auto &[edge_id, nni_info] : nni_infos) {
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": PCSPS [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << nni_info.adj_pcsps[adj_nni_type].ToHashString(5) << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": OLD_EDGES [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << !nni_info.do_optimize_edge[adj_nni_type] << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": PVS [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << GetPVs().ToHashString(nni_info.ref_primary_pv_ids[adj_nni_type],
  //       5)
  //                 << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": BEFORE_BLS [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       if (nni_info.adj_edge_ids[adj_nni_type] == NoId) continue;
  //       std::cout << DblToScientific(
  //                        branch_handler_(nni_info.adj_edge_ids[adj_nni_type]))
  //                 << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //   }
  // }

  os << "UpdateEngineAfterModifyingDAG::Leafward/RootwardPass: " << timer.Lap()
     << std::endl;

  if (IsOptimizeNewEdges()) {
    Stopwatch optimization_timer(true, Stopwatch::TimeScale::SecondScale);
    for (size_t iter = 0; iter < optimize_max_iter_; iter++) {
      // Optimize each NNI.
      for (const auto edge_id : nni_edges) {
        const auto &nni_info = nni_infos[edge_id];
        const auto &adj_edge_ids = nni_info.adj_edge_ids;
        const auto &do_optimize_edge = nni_info.do_optimize_edge;
        OptimizeEdge(NNIAdjacent::LeftChild, adj_edge_ids.left_child,
                     adj_edge_ids.focal, edge_id, false, true,
                     do_optimize_edge.left_child);
        OptimizeEdge(NNIAdjacent::RightChild, adj_edge_ids.right_child,
                     adj_edge_ids.focal, edge_id, false, true,
                     do_optimize_edge.right_child);
        OptimizeEdge(NNIAdjacent::Sister, adj_edge_ids.sister, adj_edge_ids.parent,
                     edge_id, false, true, do_optimize_edge.sister);
        OptimizeEdge(NNIAdjacent::Focal, adj_edge_ids.focal, adj_edge_ids.parent,
                     edge_id, true, true, do_optimize_edge.focal);
        if (!GetDAG().IsEdgeRoot(adj_edge_ids.parent)) {
          const auto &parent_choice = GetTPEngine().GetChoiceMap(adj_edge_ids.parent);
          OptimizeEdge(NNIAdjacent::Parent, adj_edge_ids.parent, parent_choice.parent,
                       edge_id, true, false, do_optimize_edge.parent);
        }
        os << "UpdateEngineAfterModifyingDAG::OptimizeCentralEdges: "
           << optimization_timer.Lap() << std::endl;

        // std::cout << "  === MID_MODIFYING_DAG ===" << std::endl;
        // if (is_on_watchlist) {
        //   for (auto &[edge_id, nni_info] : nni_infos) {
        //     std::cout << nni_info.post_nni.ToHashString(5) << ": PCSPS [ ";
        //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
        //       std::cout << nni_info.adj_pcsps[adj_nni_type].ToHashString(5) << " ";
        //     }
        //     std::cout << "]" << std::endl;
        //     std::cout << nni_info.post_nni.ToHashString(5) << ": PVS [ ";
        //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
        //       std::cout << GetPVs().ToHashString(
        //                        nni_info.ref_primary_pv_ids[adj_nni_type], 5)
        //                 << " ";
        //     }
        //     std::cout << "]" << std::endl;
        //     std::cout << nni_info.post_nni.ToHashString(5) << ": MID_BLS [ ";
        //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
        //       if (nni_info.adj_edge_ids[adj_nni_type] == NoId) continue;
        //       std::cout << DblToScientific(
        //                        branch_handler_(nni_info.adj_edge_ids[adj_nni_type]))
        //                 << " ";
        //     }
        //     std::cout << "]" << std::endl;
        //   }
        // }
      }
      // Optimize incidental new edges.
      for (const auto edge_id : extra_edges) {
        Stopwatch optimization_timer(true, Stopwatch::TimeScale::SecondScale);
        const auto &choice = GetTPEngine().GetChoiceMap(edge_id);
        if (!GetDAG().IsEdgeRoot(choice.parent)) {
          OptimizeEdge(NNIAdjacent::Focal, edge_id, choice.parent, edge_id);
        }
        os << "UpdateEngineAfterModifyingDAG::OptimizeAdjEdges: "
           << optimization_timer.Lap() << std::endl;
      }
      // Update new NNI PVs.
      RootwardPass();
      LeafwardPass();
      os << "UpdateEngineAfterModifyingDAG::NNIUpdatePVs: " << optimization_timer.Lap()
         << std::endl;
    }
    RootwardPass();
    LeafwardPass();
    os << "UpdateEngineAfterModifyingDAG::OptimizeNewEdges: " << timer.Lap()
       << std::endl;
  }

  // Update scores.
  // ComputeScores({{update_edges.begin(), update_edges.end()}});
  ComputeScores();
  os << "UpdateEngineAfterModifyingDAG::ComputeScores: " << timer.Lap() << std::endl;

  // std::cout << "  === AFTER_MODIFYING_DAG ===" << std::endl;
  // if (is_on_watchlist) {
  //   for (auto &[edge_id, nni_info] : nni_infos) {
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": PCSPS [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << nni_info.adj_pcsps[adj_nni_type].ToHashString(5) << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": PVS [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << GetPVs().ToHashString(nni_info.ref_primary_pv_ids[adj_nni_type],
  //       5)
  //                 << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //     std::cout << nni_info.post_nni.ToHashString(5) << ": AFTER_BLS [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       if (nni_info.adj_edge_ids[adj_nni_type] == NoId) continue;
  //       std::cout << DblToScientific(
  //                        branch_handler_(nni_info.adj_edge_ids[adj_nni_type]))
  //                 << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //   }
  // }
}

double TPEvalEngineViaLikelihood::GetTopTreeScoreWithEdge(const EdgeId edge_id) const {
  return TPEvalEngine::GetTopTreeScoreWithEdge(edge_id);
}

double TPEvalEngineViaLikelihood::GetTopTreeScoreWithProposedNNI(
    const NNIOperation &post_nni, const NNIOperation &pre_nni,
    const size_t spare_offset, std::optional<BitsetEdgeIdMap> best_edge_map_opt) {
  auto nni_info =
      GetProposedNNIInfo(post_nni, pre_nni, spare_offset, best_edge_map_opt);
  const auto &temp_pv_ids = nni_info.temp_pv_ids;
  const auto &temp_edge_ids = nni_info.temp_edge_ids;
  const auto &ref_pv_ids = nni_info.ref_pv_ids;
  const auto &ref_edge_ids = nni_info.ref_edge_ids;
  // const auto &ref_node_ids = nni_info.ref_node_ids;
  const auto &adj_edge_ids = nni_info.adj_edge_ids;
  // const auto &adj_pcsps = nni_info.adj_pcsps;
  auto &do_optimize_edge = nni_info.do_optimize_edge;

  // Initialize branch lengths.
  for (auto nni_adj_type : NNIAdjacentEnum::Iterator()) {
    // Initialize with default branch lengths.
    branch_handler_(temp_edge_ids[nni_adj_type]) =
        branch_handler_.GetDefaultBranchLength();
    // If flagged for (or best_edge_map has been given), copy branch lengths mapped over
    // from pre-edge to post-edge.
    if (do_init_proposed_branch_lengths_with_dag_ or best_edge_map_opt.has_value()) {
      branch_handler_(temp_edge_ids[nni_adj_type]) =
          branch_handler_(ref_edge_ids[nni_adj_type]);
      // If edge already exists in the DAG, use it's branch length.
      if (adj_edge_ids[nni_adj_type] != NoId) {
        branch_handler_(temp_edge_ids[nni_adj_type]) =
            branch_handler_(adj_edge_ids[nni_adj_type]);
        if (do_fix_proposed_branch_lengths_from_dag_) {
          do_optimize_edge[nni_adj_type] = false;
        }
      }
    }
  }

  // // TODO Print Branch Lengths.
  // bool is_on_watchlist = watch_all;
  // is_on_watchlist |= nni_watchlist.find(nni_info.adj_pcsps.focal.ToHashString(5)) !=
  //                    nni_watchlist.end();
  // if (is_on_watchlist) {
  //   std::cout << post_nni.ToHashString(5) << ": PCSPS [ ";
  //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //     std::cout << adj_pcsps[adj_nni_type].ToHashString(5) << " ";
  //   }
  //   std::cout << "]" << std::endl;

  //   std::cout << nni_info.post_nni.ToHashString(5) << ": ADJ_CHOICES " << std::endl;
  //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //     auto edge_id = nni_info.adj_edge_ids[adj_nni_type];
  //     if (edge_id == NoId) {
  //       std::cout << "  [ ]" << std::endl;
  //       continue;
  //     }
  //     auto pcsps = GetTPEngine().GetChoiceMap().GetEdgeChoicePCSPs(edge_id);
  //     std::cout << "  [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << pcsps[adj_nni_type].ToHashString(5) << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //   }

  //   std::cout << nni_info.post_nni.ToHashString(5) << ": REF_CHOICES " << std::endl;
  //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //     auto edge_id = nni_info.ref_edge_ids[adj_nni_type];
  //     if (edge_id == NoId) {
  //       std::cout << "  [ ]" << std::endl;
  //       continue;
  //     }
  //     auto pcsps = GetTPEngine().GetChoiceMap().GetEdgeChoicePCSPs(edge_id);
  //     std::cout << "  [ ";
  //     for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //       std::cout << pcsps[adj_nni_type].ToHashString(5) << " ";
  //     }
  //     std::cout << "]" << std::endl;
  //   }

  //   std::cout << nni_info.post_nni.ToHashString(5) << ": OLD_EDGES [ ";
  //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //     std::cout << (nni_info.adj_edge_ids[adj_nni_type] != NoId) << " ";
  //   }
  //   std::cout << "]" << std::endl;
  //   std::cout << nni_info.post_nni.ToHashString(5) << ": PVS [ ";
  //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //     std::cout << GetPVs().ToHashString(nni_info.ref_primary_pv_ids[adj_nni_type],
  //     5)
  //               << " ";
  //   }
  //   std::cout << "]" << std::endl;
  //   std::cout << post_nni.ToHashString(5) << ": BEFORE_BLS [ ";
  //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //     std::cout << DblToScientific(branch_handler_(temp_edge_ids[adj_nni_type])) << "
  //     ";
  //   }
  //   std::cout << "]" << std::endl;
  // }

  auto RootwardPass = [&]() {
    // Evolve up child P-PLVs.
    SetToEvolvedPV(temp_pv_ids.child_phatleft_, temp_edge_ids.left_child,
                   ref_pv_ids.leftchild_p_);
    SetToEvolvedPV(temp_pv_ids.child_phatright_, temp_edge_ids.right_child,
                   ref_pv_ids.rightchild_p_);
    MultiplyPVs(temp_pv_ids.child_p_, temp_pv_ids.child_phatleft_,
                temp_pv_ids.child_phatright_);
    // Evolve up parent P-PLVs
    SetToEvolvedPV(temp_pv_ids.parent_phatsister_, temp_edge_ids.sister,
                   ref_pv_ids.sister_p_);
    SetToEvolvedPV(temp_pv_ids.parent_phatfocal_, temp_edge_ids.focal,
                   temp_pv_ids.child_p_);
    MultiplyPVs(temp_pv_ids.parent_p_, temp_pv_ids.parent_phatfocal_,
                temp_pv_ids.parent_phatsister_);
  };
  auto LeafwardPass = [&]() {
    // If the parent is not the DAG root, then evolve grandparent down to parent.
    if ((ref_pv_ids.grandparent_rfocal_ != NoId)) {
      SetToEvolvedPV(temp_pv_ids.parent_rhat_, temp_edge_ids.parent,
                     ref_pv_ids.grandparent_rfocal_);
    } else {
      TakePVValue(temp_pv_ids.parent_rhat_, ref_pv_ids.parent_rhat_);
    }

    // Evolve down parent R-PLVs.
    MultiplyPVs(temp_pv_ids.parent_rfocal_, temp_pv_ids.parent_rhat_,
                temp_pv_ids.parent_phatsister_);
    MultiplyPVs(temp_pv_ids.parent_rsister_, temp_pv_ids.parent_rhat_,
                temp_pv_ids.parent_phatfocal_);
    // Evolve down child R-PLVs.
    SetToEvolvedPV(temp_pv_ids.child_rhat_, temp_edge_ids.focal,
                   temp_pv_ids.parent_rfocal_);
    MultiplyPVs(temp_pv_ids.child_rleft_, temp_pv_ids.child_rhat_,
                temp_pv_ids.child_phatright_);
    MultiplyPVs(temp_pv_ids.child_rright_, temp_pv_ids.child_rhat_,
                temp_pv_ids.child_phatleft_);
  };

  // Optimize branch lengths.
  auto OptimizeEdge = [&](size_t iter, const NNIAdjacent nni_adj_type,
                          const EdgeId edge_id, const EdgeId focal_edge_id,
                          const PVId parent_p, const PVId parent_phatfocal,
                          const PVId parent_phatsister, const PVId parent_rhat,
                          const PVId parent_rfocal, const PVId parent_rsister,
                          const PVId child_p, const PVId child_phatleft,
                          const PVId child_phatright, const bool update_branch_length,
                          const bool is_not_child_edge, const bool is_not_parent_edge) {
    if (is_not_child_edge) {
      // If child does not go to a leaf, update child_p (in case leftchild
      // or rightchild branch length were changed).
      MultiplyPVs(child_p, child_phatleft, child_phatright);
    }
    if (is_not_parent_edge) {
      // Update parent_rfocal (in case that sister branch length was
      // changed).
      MultiplyPVs(parent_rfocal, parent_rhat, parent_phatsister);
    }
    if (update_branch_length) {
      // Optimize branch length.
      // if (is_on_watchlist) {
      //   std::cout << "PROP::optimize_edge::"
      //             << nni_info.adj_pcsps[nni_adj_type].ToHashString(5) << " "
      //             << GetPVs().ToHashString(parent_rfocal, 5) << " "
      //             << GetPVs().ToHashString(child_p, 5) << std::endl;
      // }
      branch_handler_.OptimizeBranchLength(edge_id, parent_rfocal, child_p, iter > 0);
    }
    if (is_not_parent_edge) {
      // Update parent_phatfocal after changing branch length.
      SetToEvolvedPV(parent_phatfocal, edge_id, child_p);
      // If not parent edge, update parent_p after changing branch length.
      MultiplyPVs(parent_p, parent_phatfocal, parent_phatsister);
    }
  };
  auto OptimizeLeftChild = [&](const size_t iter, const bool do_optimize = true) {
    OptimizeEdge(iter, NNIAdjacent::LeftChild, temp_edge_ids.left_child,
                 temp_edge_ids.focal, temp_pv_ids.child_p_, temp_pv_ids.child_phatleft_,
                 temp_pv_ids.child_phatright_, temp_pv_ids.child_rhat_,
                 temp_pv_ids.child_rleft_, temp_pv_ids.child_rright_,
                 ref_pv_ids.leftchild_p_, PVId(NoId), PVId(NoId),
                 (do_optimize_edge.left_child && do_optimize), false, true);
  };
  auto OptimizeRightChild = [&](const size_t iter, const bool do_optimize = true) {
    OptimizeEdge(
        iter, NNIAdjacent::RightChild, temp_edge_ids.right_child, temp_edge_ids.focal,
        temp_pv_ids.child_p_, temp_pv_ids.child_phatright_, temp_pv_ids.child_phatleft_,
        temp_pv_ids.child_rhat_, temp_pv_ids.child_rright_, temp_pv_ids.child_rleft_,
        ref_pv_ids.rightchild_p_, PVId(NoId), PVId(NoId),
        (do_optimize_edge.right_child && do_optimize), false, true);
  };
  auto OptimizeSister = [&](const size_t iter, const bool do_optimize = true) {
    OptimizeEdge(iter, NNIAdjacent::Sister, temp_edge_ids.sister, temp_edge_ids.focal,
                 temp_pv_ids.parent_p_, temp_pv_ids.parent_phatsister_,
                 temp_pv_ids.parent_phatfocal_, temp_pv_ids.parent_rhat_,
                 temp_pv_ids.parent_rsister_, temp_pv_ids.parent_rfocal_,
                 ref_pv_ids.sister_p_, PVId(NoId), PVId(NoId),
                 (do_optimize_edge.sister && do_optimize), false, true);
  };
  auto OptimizeCentral = [&](const size_t iter, const bool do_optimize = true) {
    OptimizeEdge(iter, NNIAdjacent::Focal, temp_edge_ids.focal, temp_edge_ids.focal,
                 temp_pv_ids.parent_p_, temp_pv_ids.parent_phatfocal_,
                 temp_pv_ids.parent_phatsister_, temp_pv_ids.parent_rhat_,
                 temp_pv_ids.parent_rfocal_, temp_pv_ids.parent_rsister_,
                 temp_pv_ids.child_p_, temp_pv_ids.child_phatleft_,
                 temp_pv_ids.child_phatright_, (do_optimize_edge.focal && do_optimize),
                 true, true);
  };
  auto OptimizeParent = [&](const size_t iter, const bool do_optimize = true) {
    OptimizeEdge(iter, NNIAdjacent::Parent, temp_edge_ids.parent, temp_edge_ids.focal,
                 PVId(NoId), PVId(NoId), PVId(NoId), ref_pv_ids.grandparent_rhat_,
                 ref_pv_ids.grandparent_rfocal_, ref_pv_ids.grandparent_rsister_,
                 temp_pv_ids.parent_p_, temp_pv_ids.parent_phatfocal_,
                 temp_pv_ids.parent_phatsister_,
                 (do_optimize_edge.parent && do_optimize), true, false);
  };

  // Initialize new PLVs.
  RootwardPass();
  LeafwardPass();

  // Branch length optimization.
  if (IsOptimizeNewEdges()) {
    for (size_t iter = 0; iter < optimize_max_iter_; iter++) {
      // Optimize Branch Lengths
      OptimizeLeftChild(iter);
      OptimizeRightChild(iter);
      OptimizeSister(iter);
      OptimizeCentral(iter);
      if (!post_nni.GetParent().SubsplitIsRootsplit()) {
        if (temp_edge_ids.parent != NoId) {
          OptimizeParent(iter);
        }
      }

      // // TODO Print NNI Info
      // if (is_on_watchlist) {
      //   std::cout << nni_info.post_nni.ToHashString(5) << ": MID_PVS [ ";
      //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
      //     std::cout <<
      //     GetPVs().ToHashString(nni_info.ref_primary_pv_ids[adj_nni_type],
      //                                        5)
      //               << " ";
      //   }
      //   std::cout << "]" << std::endl;
      //   std::cout << post_nni.ToHashString(5) << ": MID_BLS [ ";
      //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
      //     std::cout << DblToScientific(branch_handler_(temp_edge_ids[adj_nni_type]))
      //               << " ";
      //   }
      //   std::cout << "]" << std::endl;
      // }

      // Update PLVs.
      RootwardPass();
      LeafwardPass();
    }
  }

  // // TODO Print NNI Info
  // if (is_on_watchlist) {
  //   std::cout << nni_info.post_nni.ToHashString(5) << ": AFTER_PVS [ ";
  //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //     std::cout << GetPVs().ToHashString(nni_info.ref_primary_pv_ids[adj_nni_type],
  //     5)
  //               << " ";
  //   }
  //   std::cout << post_nni.ToHashString(5) << ": AFTER_BLS [ ";
  //   for (auto adj_nni_type : NNIAdjacentEnum::Iterator()) {
  //     std::cout << DblToScientific(branch_handler_(temp_edge_ids[adj_nni_type])) << "
  //     ";
  //   }
  //   std::cout << "]" << std::endl;
  // }

  // Compute likelihood of focal edge by evolving up from child to parent.
  ComputeLikelihood(temp_edge_ids.focal, temp_pv_ids.child_p_,
                    temp_pv_ids.parent_rfocal_);
  EigenVectorXd top_tree_likelihood =
      GetMatrix().block(temp_edge_ids.focal.value_, 0, 1, GetMatrix().cols()) *
      GetTPEngine().GetSitePatternWeights();
  return top_tree_likelihood[0];
}

ProposedNNIInfo TPEvalEngineViaLikelihood::GetProposedNNIInfo(
    const NNIOperation &post_nni, const NNIOperation &tmp_pre_nni,
    const size_t spare_offset, std::optional<BitsetEdgeIdMap> best_edge_map_opt) const {
  // Get the best pre-NNI candidate.
  const NNIOperation pre_nni =
      GetTPEngine().FindHighestPriorityNeighborNNIInDAG(post_nni);
  // Get reference edge choices from pre-NNI in DAG, remapped according to post-NNI.
  const auto pre_edge_id = GetDAG().GetEdgeIdx(pre_nni);
  const auto clade_map =
      NNIOperation::BuildNNICladeMapFromPreNNIToNNI(pre_nni, post_nni);
  const auto pre_edge_ids = GetTPEngine().GetChoiceMap().GetEdgeChoice(pre_edge_id);
  const auto pre_edge_ids_remapped =
      GetTPEngine().GetChoiceMap().RemapEdgeChoiceDataViaNNICladeMap(pre_edge_ids,
                                                                     clade_map);
  NNIAdjEdgeIds ref_edge_ids;
  ref_edge_ids.parent = pre_edge_ids_remapped.parent;
  ref_edge_ids.sister = pre_edge_ids_remapped.sister;
  ref_edge_ids.focal = pre_edge_id;
  ref_edge_ids.left_child = pre_edge_ids_remapped.left_child;
  ref_edge_ids.right_child = pre_edge_ids_remapped.right_child;
  // Get reference node choices from pre-NNI in DAG, remapped according to post-NNI.
  const auto pre_node_ids_remapped =
      GetTPEngine().GetChoiceMap().GetEdgeChoiceNodeIds(pre_edge_ids_remapped);
  NNIAdjNodeIds ref_node_ids;
  ref_node_ids.parent = pre_node_ids_remapped.parent;
  ref_node_ids.sister = pre_node_ids_remapped.sister;
  ref_node_ids.focal = NodeId(NoId);
  ref_node_ids.left_child = pre_node_ids_remapped.left_child;
  ref_node_ids.right_child = pre_node_ids_remapped.right_child;

  // Get adjacent PCSPs from edges.
  const auto adj_pcsps =
      GetTPEngine().BuildAdjacentPCSPsToProposedNNI(post_nni, pre_node_ids_remapped);
  // If passed best_edge_map, override ref_edge_ids using map.
  if (best_edge_map_opt.has_value()) {
    auto &best_edge_map = best_edge_map_opt.value();
    for (auto nni_adj_type : NNIAdjacentEnum::Iterator()) {
      ref_edge_ids[nni_adj_type] = best_edge_map[adj_pcsps[nni_adj_type]];
    }
  }
  // Update adjacent edge ids if they already exists in DAG.
  NNIAdjEdgeIds adj_edge_ids;
  for (auto nni_adj_type : NNIAdjacentEnum::Iterator()) {
    EdgeId edge_id{NoId};
    if (GetDAG().ContainsEdge(adj_pcsps[nni_adj_type])) {
      edge_id = GetDAG().GetEdgeIdx(adj_pcsps[nni_adj_type]);
    }
    adj_edge_ids[nni_adj_type] = edge_id;
  }
  NNIAdjBools do_optimize_edge{true, true, true, true, true};

  // Get reference PLV IDs from pre-NNI in DAG, remapped according to post-NNI.
  const auto pre_pv_ids = GetLocalPVIdsOfEdge(pre_edge_id);
  const auto pre_pv_ids_remapped = RemapLocalPVIdsForPostNNI(pre_pv_ids, clade_map);
  NNIAdjacentMap<PVId> ref_primary_pv_ids;
  ref_primary_pv_ids.parent = pre_pv_ids_remapped.grandparent_rfocal_;
  ref_primary_pv_ids.sister = pre_pv_ids_remapped.sister_p_;
  ref_primary_pv_ids.focal = pre_pv_ids_remapped.parent_rfocal_;
  ref_primary_pv_ids.left_child = pre_pv_ids_remapped.leftchild_p_;
  ref_primary_pv_ids.right_child = pre_pv_ids_remapped.rightchild_p_;
  // Get temp locations for post-NNI PVs.
  const auto temp_pv_ids = GetTempLocalPVIdsForProposedNNIs(spare_offset);
  // Get temp locations for post-NNI branch lengths.
  const auto temp_edge_ids = GetTempEdgeIdsForProposedNNIs(spare_offset);

  ProposedNNIInfo nni_info;
  nni_info.post_nni = post_nni;
  nni_info.pre_nni = pre_nni;
  nni_info.temp_pv_ids = temp_pv_ids;
  nni_info.temp_edge_ids = temp_edge_ids;
  nni_info.ref_pv_ids = pre_pv_ids_remapped;
  nni_info.ref_primary_pv_ids = ref_primary_pv_ids;
  nni_info.ref_edge_ids = ref_edge_ids;
  nni_info.ref_node_ids = ref_node_ids;
  nni_info.adj_edge_ids = adj_edge_ids;
  nni_info.adj_pcsps = adj_pcsps;
  nni_info.do_optimize_edge = do_optimize_edge;

  return nni_info;
}

ProposedNNIInfo TPEvalEngineViaLikelihood::GetRealNNIInfo(
    const NNIOperation &post_nni, const NNIOperation &tmp_pre_nni,
    const size_t spare_offset, std::optional<BitsetEdgeIdMap> best_edge_map_opt) const {
  auto nni_info =
      GetProposedNNIInfo(post_nni, tmp_pre_nni, spare_offset, best_edge_map_opt);

  auto edge_id = GetDAG().GetEdgeIdx(post_nni);
  nni_info.ref_pv_ids = GetLocalPVIdsOfEdge(edge_id);
  auto &ref_pv_ids = nni_info.ref_pv_ids;
  auto &ref_primary_pv_ids = nni_info.ref_primary_pv_ids;
  ref_primary_pv_ids.parent = ref_pv_ids.grandparent_rfocal_;
  ref_primary_pv_ids.sister = ref_pv_ids.sister_p_;
  ref_primary_pv_ids.focal = ref_pv_ids.parent_rfocal_;
  ref_primary_pv_ids.left_child = ref_pv_ids.leftchild_p_;
  ref_primary_pv_ids.right_child = ref_pv_ids.rightchild_p_;

  return nni_info;
}

LocalPVIds TPEvalEngineViaLikelihood::GetTempLocalPVIdsForProposedNNIs(
    const size_t spare_offset) const {
  LocalPVIds temp_pv_ids;
  size_t spare_count = 0;
  auto GetNextSparePVIndex = [this, spare_offset, &spare_count]() {
    PVId next_pvid =
        GetPVs().GetSparePVIndex((spare_offset * spare_nodes_per_nni_) + spare_count);
    spare_count++;
    return next_pvid;
  };

  temp_pv_ids.grandparent_rhat_ = GetNextSparePVIndex();
  temp_pv_ids.grandparent_rfocal_ = GetNextSparePVIndex();
  temp_pv_ids.grandparent_rsister_ = GetNextSparePVIndex();

  temp_pv_ids.parent_p_ = GetNextSparePVIndex();
  temp_pv_ids.parent_phatfocal_ = GetNextSparePVIndex();
  temp_pv_ids.parent_phatsister_ = GetNextSparePVIndex();
  temp_pv_ids.parent_rfocal_ = GetNextSparePVIndex();
  temp_pv_ids.parent_rhat_ = GetNextSparePVIndex();
  temp_pv_ids.parent_rsister_ = GetNextSparePVIndex();

  temp_pv_ids.child_p_ = GetNextSparePVIndex();
  temp_pv_ids.child_phatleft_ = GetNextSparePVIndex();
  temp_pv_ids.child_phatright_ = GetNextSparePVIndex();
  temp_pv_ids.child_rhat_ = GetNextSparePVIndex();
  temp_pv_ids.child_rleft_ = GetNextSparePVIndex();
  temp_pv_ids.child_rright_ = GetNextSparePVIndex();

  temp_pv_ids.sister_p_ = GetNextSparePVIndex();
  temp_pv_ids.leftchild_p_ = GetNextSparePVIndex();
  temp_pv_ids.rightchild_p_ = GetNextSparePVIndex();

  return temp_pv_ids;
}

NNIAdjEdgeIds TPEvalEngineViaLikelihood::GetTempEdgeIdsForProposedNNIs(
    const size_t spare_offset) const {
  NNIAdjEdgeIds temp_edge_ids;
  size_t spare_count = 0;
  auto GetNextSpareEdgeId = [this, spare_offset, &spare_count]() {
    EdgeId next_edge_id = EdgeId(GetTPEngine().GetEdgeCount() +
                                 (spare_offset * spare_edges_per_nni_) + spare_count);
    spare_count++;
    return next_edge_id;
  };
  temp_edge_ids.focal = GetNextSpareEdgeId();
  temp_edge_ids.parent = GetNextSpareEdgeId();
  temp_edge_ids.sister = GetNextSpareEdgeId();
  temp_edge_ids.left_child = GetNextSpareEdgeId();
  temp_edge_ids.right_child = GetNextSpareEdgeId();
  return temp_edge_ids;
}

void TPEvalEngineViaLikelihood::CopyEdgeData(const EdgeId src_edge_id,
                                             const EdgeId dest_edge_id) {
  branch_handler_(dest_edge_id) = branch_handler_(src_edge_id);
  TPEvalEngine::CopyEdgeData(src_edge_id, dest_edge_id);
}

// For rootward traversal. Compute the likelihood PV for a given node, by using
// the left and right child edges from the choice map of the edge below node.
void TPEvalEngineViaLikelihood::PopulateRootwardPVForNode(const NodeId node_id) {
  for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto adj_node_id :
         GetDAG().GetDAGNode(node_id).GetNeighbors(Direction::Rootward, clade)) {
      const EdgeId edge_id = GetDAG().GetEdgeIdx(adj_node_id, node_id);
      PopulateRootwardPVForEdge(edge_id);
    }
  }
}

void TPEvalEngineViaLikelihood::PopulateRootwardPVForEdge(const EdgeId edge_id) {
  // Populate edge PLV by evolving up given edge.
  const auto choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(edge_id);
  // Update parent's PLeftHat PLV, evolved up from left child's P PLV.
  if (choices.left_child != NoId) {
    EvolvePPVUpEdge(edge_id, choices.left_child);
  }
  // Update parent's PRightHat PLV, evolved up from right child's P PLV.
  if (choices.right_child != NoId) {
    EvolvePPVUpEdge(edge_id, choices.right_child);
  }
  // Update P-PLV from PHatLeft and PHatRight by taking product.
  const PVId p_pvid = GetPVs().GetPVIndex(PLVType::P, edge_id);
  const PVId phatleft_pvid = GetPVs().GetPVIndex(PLVType::PHatLeft, edge_id);
  const PVId phatright_pvid = GetPVs().GetPVIndex(PLVType::PHatRight, edge_id);
  if ((choices.left_child != NoId) && (choices.right_child != NoId)) {
    MultiplyPVs(p_pvid, phatleft_pvid, phatright_pvid);
  } else if (choices.left_child != NoId) {
    TakePVValue(p_pvid, phatleft_pvid);
  } else if (choices.right_child != NoId) {
    TakePVValue(p_pvid, phatright_pvid);
  }
}

// For leafward traversal. Compute the likelihood PV for a given node, by using the
// parent and sister edges from the choice map of the edge above node.
void TPEvalEngineViaLikelihood::PopulateLeafwardPVForNode(const NodeId node_id) {
  for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto adj_node_id :
         GetDAG().GetDAGNode(node_id).GetNeighbors(Direction::Leafward, clade)) {
      const EdgeId edge_id = GetDAG().GetEdgeIdx(node_id, adj_node_id);
      PopulateLeafwardPVForEdge(edge_id);
    }
  }
}

void TPEvalEngineViaLikelihood::PopulateLeafwardPVForEdge(const EdgeId edge_id) {
  // Populate edge PLV by evolving down given edge.
  const auto choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(edge_id);
  // Evolve down parent.
  // Updates child's R-PLV by evolving down from parent's RFocalHat-PLV.
  if (choices.parent != NoId) {
    EvolveRPVDownEdge(choices.parent, edge_id);
  }
  // Updates child's RFocalHat-PLV and RSisterHat-PLV by taking (RHat o PHatOpposite).
  const PVId rhat_pvid = GetPVs().GetPVIndex(PLVType::RHat, edge_id);
  const PVId rleft_pvid = GetPVs().GetPVIndex(PLVType::RLeft, edge_id);
  const PVId rright_pvid = GetPVs().GetPVIndex(PLVType::RRight, edge_id);
  const PVId phatleft_pvid = GetPVs().GetPVIndex(PLVType::PHatLeft, edge_id);
  const PVId phatright_pvid = GetPVs().GetPVIndex(PLVType::PHatRight, edge_id);
  MultiplyPVs(rleft_pvid, rhat_pvid, phatright_pvid);
  MultiplyPVs(rright_pvid, rhat_pvid, phatleft_pvid);
}

void TPEvalEngineViaLikelihood::PopulateLeafPVsWithSitePatterns(
    std::optional<EdgeIdVector> opt_edge_ids) {
  // TODO add option to give specific edge_ids.
  std::ignore = opt_edge_ids;

  auto BuildPVForSitePattern = [this](const TaxonId taxon_id, const PVId pv_id) {
    auto &pv = GetPVs().GetPV(pv_id);
    pv.setZero();
    const auto &pattern = GetSitePattern().GetPatterns()[taxon_id.value_];
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative symbol!");
      if (symbol == MmappedNucleotidePLV::base_count_) {  // Gap character.
        pv.col(site_idx).setConstant(1.);
      } else if (symbol < MmappedNucleotidePLV::base_count_) {
        pv(symbol, site_idx) = 1.;
      }
      site_idx++;
    }
  };

  for (const auto taxon_id : GetDAG().GetTaxonIds()) {
    const EdgeIdVector leaf_edge_ids = GetDAG().GetLeafEdgeIds(taxon_id);
    PVId src_pvid;
    if (leaf_edge_ids.size() > 0) {
      src_pvid = GetPVs().GetPVIndex(PLVType::P, leaf_edge_ids[0]);
      BuildPVForSitePattern(taxon_id, src_pvid);
    }
    for (const auto edge_id : leaf_edge_ids) {
      for (const auto pv_type : {PLVType::P}) {
        const auto dest_pvid = GetPVs().GetPVIndex(pv_type, edge_id);
        TakePVValue(dest_pvid, src_pvid);
      }
    }
  }
}

void TPEvalEngineViaLikelihood::PopulateRootPVsWithStationaryDistribution(
    std::optional<EdgeIdVector> opt_edge_ids) {
  auto AssignPVToStationaryDistribution = [this](const PVId pv_id) {
    auto &pv = GetPVs().GetPV(pv_id);
    for (Eigen::Index row_idx = 0; row_idx < pv.rows(); ++row_idx) {
      pv.row(row_idx).array() = stationary_distribution_(row_idx);
    }
  };

  EdgeIdVector rootsplit_edge_ids =
      (opt_edge_ids.has_value() ? opt_edge_ids.value()
                                : GetDAG().GetRootsplitEdgeIds());
  for (const auto edge_id : rootsplit_edge_ids) {
    for (const auto pv_type : {PLVType::RHat}) {
      const auto pvid = GetPVs().GetPVIndex(pv_type, edge_id);
      AssignPVToStationaryDistribution(pvid);
    }
  }
}

void TPEvalEngineViaLikelihood::ComputeScores(
    std::optional<EdgeIdVector> opt_edge_ids) {
  EdgeIdVector edge_ids = opt_edge_ids.has_value()
                              ? opt_edge_ids.value()
                              : GetDAG().LeafwardEdgeTraversalTrace(true);
  // Compute likelihoods for each edge.
  for (const auto edge_id : edge_ids) {
    // const auto choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(edge_id);
    const auto &[parent_pvid, child_pvid] = GetPrimaryPVIdsOfEdge(edge_id);
    ComputeLikelihood(edge_id, child_pvid, parent_pvid);
  }
  auto &top_tree_likelihoods = GetTopTreeScores();
  top_tree_likelihoods =
      GetMatrix().block(0, 0, GetTPEngine().GetEdgeCount(), GetMatrix().cols()) *
      GetTPEngine().GetSitePatternWeights();
}

void TPEvalEngineViaLikelihood::InitializeBranchLengthHandler() {
  // Set Nongradient Brent.
  DAGBranchHandler::NegLogLikelihoodFunc brent_nongrad_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double log_branch_length) {
        SetTransitionMatrixToHaveBranchLength(exp(log_branch_length));
        PreparePerPatternLogLikelihoodsForEdge(parent_id, child_id);
        double result =
            -per_pattern_log_likelihoods_.dot(GetTPEngine().GetSitePatternWeights());
        return result;
      };
  branch_handler_.SetBrentFunc(brent_nongrad_func);
  // Set Gradient Brent.
  DAGBranchHandler::NegLogLikelihoodAndDerivativeFunc brent_grad_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double log_branch_length) {
        double branch_length = exp(log_branch_length);
        branch_handler_(edge_id) = branch_length;
        auto [log_likelihood, log_likelihood_derivative] =
            this->LogLikelihoodAndDerivative(edge_id);
        return std::make_pair(-log_likelihood,
                              -branch_length * log_likelihood_derivative);
      };
  branch_handler_.SetBrentWithGradientFunc(brent_grad_func);
  // Set Gradient Ascent.
  DAGBranchHandler::LogLikelihoodAndDerivativeFunc grad_ascent_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double branch_length) {
        branch_handler_(edge_id) = branch_length;
        return this->LogLikelihoodAndDerivative(edge_id);
      };
  branch_handler_.SetGradientAscentFunc(grad_ascent_func);
  // Set Logspace Gradient Ascent.
  DAGBranchHandler::LogLikelihoodAndDerivativeFunc logspace_grad_ascent_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double branch_length) {
        branch_handler_(edge_id) = branch_length;
        return this->LogLikelihoodAndDerivative(edge_id);
      };
  branch_handler_.SetLogSpaceGradientAscentFunc(logspace_grad_ascent_func);
  // Set Newton-Raphson.
  DAGBranchHandler::LogLikelihoodAndFirstTwoDerivativesFunc newton_raphson_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double log_branch_length) {
        double x = exp(log_branch_length);
        branch_handler_(edge_id) = x;
        auto [f_x, f_prime_x, f_double_prime_x] =
            this->LogLikelihoodAndFirstTwoDerivatives(edge_id);
        // x = exp(y) --> f'(exp(y)) = exp(y) * f'(exp(y)) = x * f'(x)
        double f_prime_y = x * f_prime_x;
        double f_double_prime_y = f_prime_y + std::pow(x, 2) * f_double_prime_x;
        return std::make_tuple(f_x, f_prime_y, f_double_prime_y);
      };
  branch_handler_.SetNewtonRaphsonFunc(newton_raphson_func);
}

void TPEvalEngineViaLikelihood::BranchLengthOptimization(
    std::optional<bool> check_branch_convergence) {
  bool check_branch_convergence_ = check_branch_convergence.has_value()
                                       ? check_branch_convergence.value()
                                       : !IsFirstOptimization();
  // Update R-PVs and optimize branch lengths leafward.
  const EdgeIdVector edge_ids = GetDAG().RootwardEdgeTraversalTrace(false);
  for (size_t opt_count = 0; opt_count < GetOptimizationMaxIteration(); opt_count++) {
    for (const EdgeId edge_id : edge_ids) {
      BranchLengthOptimization(edge_id, check_branch_convergence_);
    }
    IncrementOptimizationCount();
  }
}

void TPEvalEngineViaLikelihood::BranchLengthOptimization(
    const EdgeId edge_id, const bool check_branch_convergence, const bool update_only) {
  const auto &choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(edge_id);
  if (choices.parent != NoId) {
    // Update parent's RFocal PLV from previous branch length changes.
    PopulateRootwardPVForEdge(edge_id);
    PopulateRootwardPVForEdge(choices.parent);
    PopulateLeafwardPVForEdge(choices.parent);
    // Optimize branch length.
    if (!update_only) {
      const auto &[parent_rfocal_pvid, child_p_pvid] = GetPrimaryPVIdsOfEdge(edge_id);
      if (parent_rfocal_pvid != NoId) {
        branch_handler_.OptimizeBranchLength(edge_id, parent_rfocal_pvid, child_p_pvid,
                                             check_branch_convergence);
      }
    }
    // Update parent's PHatFocal PLV for new branch length.
    PopulateLeafwardPVForEdge(edge_id);
  }
}

void TPEvalEngineViaLikelihood::EvolvePPVUpEdge(const EdgeId rootward_edge_id,
                                                const EdgeId leafward_edge_id) {
  const auto focal = GetDAG().GetFocalClade(leafward_edge_id);
  const PVId parent_phatfocal_pvid =
      GetPVs().GetPVIndex(PLVTypeEnum::PPLVType(focal), rootward_edge_id);
  const PVId child_p_pvid = GetPVs().GetPVIndex(PLVType::P, leafward_edge_id);
  SetToEvolvedPV(parent_phatfocal_pvid, leafward_edge_id, child_p_pvid);
}

void TPEvalEngineViaLikelihood::EvolveRPVDownEdge(const EdgeId rootward_edge_id,
                                                  const EdgeId leafward_edge_id) {
  const auto focal = GetDAG().GetFocalClade(leafward_edge_id);
  const PVId parent_rfocal_pvid =
      GetPVs().GetPVIndex(PLVTypeEnum::RPLVType(focal), rootward_edge_id);
  const PVId child_rhat_pvid = GetPVs().GetPVIndex(PLVType::RHat, leafward_edge_id);
  SetToEvolvedPV(child_rhat_pvid, leafward_edge_id, parent_rfocal_pvid);
}

std::pair<PVId, PVId> TPEvalEngineViaLikelihood::GetPrimaryPVIdsOfEdge(
    const EdgeId edge_id) const {
  const auto &choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(edge_id);
  PVId parent_rfocal_pvid{NoId};
  if (choices.parent == NoId) {
    // return std::make_pair(PVId(NoId), PVId(NoId));
    EdgeId root_id = GetDAG().GetFirstRootsplitEdgeId();
    parent_rfocal_pvid = GetPVs().GetPVIndex(PLVType::RHat, root_id);
  } else {
    parent_rfocal_pvid = GetPVs().GetPVIndex(
        PLVTypeEnum::RPLVType(GetDAG().GetFocalClade(edge_id)), choices.parent);
  }
  PVId child_p_pvid = GetPVs().GetPVIndex(PLVType::P, edge_id);
  return std::make_pair(parent_rfocal_pvid, child_p_pvid);
}

LocalPVIds TPEvalEngineViaLikelihood::GetLocalPVIdsOfEdge(const EdgeId edge_id) const {
  LocalPVIds pv_ids;
  const auto &choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(edge_id);
  // Grandparent PVIds (if parent edge is not a rootsplit).
  pv_ids.grandparent_rhat_ = GetPVs().GetPVIndex(PLVType::RHat, choices.parent);
  if (!GetDAG().IsEdgeRoot(choices.parent)) {
    const auto &parent_choices =
        GetTPEngine().GetChoiceMap().GetEdgeChoice(choices.parent);
    const auto focal = GetDAG().GetFocalClade(choices.parent);
    const auto sister = Bitset::Opposite(focal);
    // pv_ids.grandparent_rhat_ =
    //     GetPVs().GetPVIndex(PLVType::RHat, parent_choices.parent);
    pv_ids.grandparent_rfocal_ =
        GetPVs().GetPVIndex(PLVTypeEnum::RPLVType(focal), parent_choices.parent);
    pv_ids.grandparent_rsister_ =
        GetPVs().GetPVIndex(PLVTypeEnum::RPLVType(sister), parent_choices.parent);
  }

  // Parent PVIds
  pv_ids.parent_p_ = GetPVs().GetPVIndex(PLVType::P, choices.parent);
  pv_ids.parent_phatfocal_ = GetPVs().GetPVIndex(
      PLVTypeEnum::PPLVType(GetDAG().GetFocalClade(edge_id)), choices.parent);
  pv_ids.parent_phatsister_ = GetPVs().GetPVIndex(
      PLVTypeEnum::PPLVType(GetDAG().GetSisterClade(edge_id)), choices.parent);
  pv_ids.parent_rhat_ = GetPVs().GetPVIndex(PLVType::RHat, choices.parent);
  pv_ids.parent_rfocal_ = GetPVs().GetPVIndex(
      PLVTypeEnum::RPLVType(GetDAG().GetFocalClade(edge_id)), choices.parent);
  pv_ids.parent_rsister_ = GetPVs().GetPVIndex(
      PLVTypeEnum::RPLVType(GetDAG().GetSisterClade(edge_id)), choices.parent);
  // Child PVIds
  pv_ids.child_p_ = GetPVs().GetPVIndex(PLVType::P, edge_id);
  pv_ids.child_phatleft_ = GetPVs().GetPVIndex(PLVType::PHatLeft, edge_id);
  pv_ids.child_phatright_ = GetPVs().GetPVIndex(PLVType::PHatRight, edge_id);
  pv_ids.child_rhat_ = GetPVs().GetPVIndex(PLVType::RHat, edge_id);
  pv_ids.child_rleft_ = GetPVs().GetPVIndex(PLVType::RLeft, edge_id);
  pv_ids.child_rright_ = GetPVs().GetPVIndex(PLVType::RRight, edge_id);
  // Grandchild PVIds
  pv_ids.sister_p_ = GetPVs().GetPVIndex(PLVType::P, choices.sister);
  pv_ids.leftchild_p_ = GetPVs().GetPVIndex(PLVType::P, choices.left_child);
  pv_ids.rightchild_p_ = GetPVs().GetPVIndex(PLVType::P, choices.right_child);

  return pv_ids;
}

LocalPVIds TPEvalEngineViaLikelihood::RemapLocalPVIdsForPostNNI(
    const LocalPVIds &pre_pv_ids, const NNIOperation::NNICladeArray &clade_map) const {
  using NNIClade = NNIOperation::NNIClade;
  using NNICladeEnum = NNIOperation::NNICladeEnum;
  LocalPVIds post_pv_ids(pre_pv_ids);
  NNICladeEnum::Array<PVId> pre_id_map;

  pre_id_map[clade_map[NNIClade::ParentSister]] = pre_pv_ids.sister_p_;
  pre_id_map[clade_map[NNIClade::ChildLeft]] = pre_pv_ids.leftchild_p_;
  pre_id_map[clade_map[NNIClade::ChildRight]] = pre_pv_ids.rightchild_p_;
  post_pv_ids.sister_p_ = pre_id_map[NNIClade::ParentSister];
  post_pv_ids.leftchild_p_ = pre_id_map[NNIClade::ChildLeft];
  post_pv_ids.rightchild_p_ = pre_id_map[NNIClade::ChildRight];

  pre_id_map[clade_map[NNIClade::ParentSister]] = pre_pv_ids.parent_rsister_;
  pre_id_map[clade_map[NNIClade::ChildLeft]] = pre_pv_ids.child_rleft_;
  pre_id_map[clade_map[NNIClade::ChildRight]] = pre_pv_ids.child_rright_;
  post_pv_ids.parent_rsister_ = pre_id_map[NNIClade::ParentSister];
  post_pv_ids.child_rleft_ = pre_id_map[NNIClade::ChildLeft];
  post_pv_ids.child_rright_ = pre_id_map[NNIClade::ChildRight];

  return post_pv_ids;
}

void TPEvalEngineViaLikelihood::TakePVValue(const PVId dest_id, const PVId src_id) {
  GetPVs().GetPV(dest_id).array() = GetPVs().GetPV(src_id).array();
}

void TPEvalEngineViaLikelihood::MultiplyPVs(const PVId dest_id, const PVId src1_id,
                                            const PVId src2_id) {
  GetPVs().GetPV(dest_id).array() =
      GetPVs().GetPV(src1_id).array() * GetPVs().GetPV(src2_id).array();
  // #462: Need to add rescaling to PVs.
}

void TPEvalEngineViaLikelihood::ComputeLikelihood(const EdgeId dest_id,
                                                  const PVId child_id,
                                                  const PVId parent_id) {
  SetTransitionMatrixToHaveBranchLength(branch_handler_(dest_id));
  PreparePerPatternLogLikelihoodsForEdge(parent_id, child_id);
  log_likelihoods_.row(dest_id.value_) = per_pattern_log_likelihoods_;
}

void TPEvalEngineViaLikelihood::SetToEvolvedPV(const PVId dest_id, const EdgeId edge_id,
                                               const PVId src_id) {
  SetTransitionMatrixToHaveBranchLength(branch_handler_(edge_id));
  GetPVs().GetPV(dest_id).array() = (transition_matrix_ * GetPVs().GetPV(src_id));
}

void TPEvalEngineViaLikelihood::MultiplyWithEvolvedPV(const PVId dest_id,
                                                      const EdgeId edge_id,
                                                      const PVId src_id) {
  SetTransitionMatrixToHaveBranchLength(branch_handler_(edge_id));
  GetPVs().GetPV(dest_id).array() =
      GetPVs().GetPV(dest_id).array() *
      (transition_matrix_ * GetPVs().GetPV(src_id)).array();
}

DoublePair TPEvalEngineViaLikelihood::LogLikelihoodAndDerivative(const EdgeId edge_id) {
  const auto &[parent_pvid, child_pvid] = GetPrimaryPVIdsOfEdge(edge_id);
  SetTransitionAndDerivativeMatricesToHaveBranchLength(branch_handler_(edge_id));
  PreparePerPatternLogLikelihoodsForEdge(parent_pvid, child_pvid);
  // The prior is expressed using the current value of q_.
  // The phylogenetic component of the likelihood is weighted with the number of times
  // we see the site patterns.
  const double log_likelihood =
      per_pattern_log_likelihoods_.dot(GetTPEngine().GetSitePatternWeights());

  // The per-site likelihood derivative is calculated in the same way as the per-site
  // likelihood, but using the derivative matrix instead of the transition matrix.
  // We first prepare two useful vectors _without_ likelihood rescaling, because the
  // rescalings cancel out in the ratio below.
  PrepareUnrescaledPerPatternLikelihoodDerivatives(parent_pvid, child_pvid);
  PrepareUnrescaledPerPatternLikelihoods(parent_pvid, child_pvid);
  // If l_i is the per-site likelihood, the derivative of log(l_i) is the derivative
  // of l_i divided by l_i.
  per_pattern_likelihood_derivative_ratios_ =
      per_pattern_likelihood_derivatives_.array() / per_pattern_likelihoods_.array();
  const double log_likelihood_derivative =
      per_pattern_likelihood_derivative_ratios_.dot(
          GetTPEngine().GetSitePatternWeights());
  return {log_likelihood, log_likelihood_derivative};
}

std::tuple<double, double, double>
TPEvalEngineViaLikelihood::LogLikelihoodAndFirstTwoDerivatives(const EdgeId edge_id) {
  const auto &[parent_pvid, child_pvid] = GetPrimaryPVIdsOfEdge(edge_id);
  SetTransitionAndDerivativeMatricesToHaveBranchLength(branch_handler_(edge_id));
  PreparePerPatternLogLikelihoodsForEdge(parent_pvid, child_pvid);

  const double log_likelihood =
      per_pattern_log_likelihoods_.dot(GetTPEngine().GetSitePatternWeights());

  // The per-site likelihood derivative is calculated in the same way as the per-site
  // likelihood, but using the derivative matrix instead of the transition matrix.
  // We first prepare two useful vectors _without_ likelihood rescaling, because the
  // rescalings cancel out in the ratio below.
  PrepareUnrescaledPerPatternLikelihoodDerivatives(parent_pvid, child_pvid);
  PrepareUnrescaledPerPatternLikelihoods(parent_pvid, child_pvid);
  // If l_i is the per-site likelihood, the derivative of log(l_i) is the derivative
  // of l_i divided by l_i.
  per_pattern_likelihood_derivative_ratios_ =
      per_pattern_likelihood_derivatives_.array() / per_pattern_likelihoods_.array();
  const double log_likelihood_gradient = per_pattern_likelihood_derivative_ratios_.dot(
      GetTPEngine().GetSitePatternWeights());
  // Second derivative is calculated the same way, but has an extra term due to
  // the product rule.
  PrepareUnrescaledPerPatternLikelihoodSecondDerivatives(parent_pvid, child_pvid);

  per_pattern_likelihood_second_derivative_ratios_ =
      (per_pattern_likelihood_second_derivatives_.array() *
           per_pattern_likelihoods_.array() -
       per_pattern_likelihood_derivatives_.array() *
           per_pattern_likelihood_derivatives_.array()) /
      (per_pattern_likelihoods_.array() * per_pattern_likelihoods_.array());

  const double log_likelihood_hessian =
      per_pattern_likelihood_second_derivative_ratios_.dot(
          GetTPEngine().GetSitePatternWeights());

  return std::make_tuple(log_likelihood, log_likelihood_gradient,
                         log_likelihood_hessian);
}

void TPEvalEngineViaLikelihood::SetTransitionMatrixToHaveBranchLength(
    double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void TPEvalEngineViaLikelihood::SetTransitionAndDerivativeMatricesToHaveBranchLength(
    double branch_length) {
  diagonal_vector_ = (branch_length * eigenvalues_).array().exp();
  diagonal_matrix_.diagonal() = diagonal_vector_;
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
  // Now calculating derivative matrix
  diagonal_matrix_.diagonal() = eigenvalues_.array() * diagonal_vector_.array();
  derivative_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
  // Now calculating hessian matrix
  diagonal_matrix_.diagonal() =
      eigenvalues_.array() * eigenvalues_.array() * diagonal_vector_.array();
  hessian_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void TPEvalEngineViaLikelihood::SetTransitionMatrixToHaveBranchLengthAndTranspose(
    double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ =
      inverse_eigenmatrix_.transpose() * diagonal_matrix_ * eigenmatrix_.transpose();
}

// ** TPEvalEngineViaParsimony

TPEvalEngineViaParsimony::TPEvalEngineViaParsimony(TPEngine &tp_engine,
                                                   const std::string &mmap_path)
    : TPEvalEngine(tp_engine),
      parsimony_pvs_(mmap_path, GetDAG().EdgeCountWithLeafSubsplits(),
                     GetSitePattern().PatternCount(), 2.0) {
  GrowNodeData(GetDAG().NodeCount(), std::nullopt, std::nullopt, true);
  GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits(), std::nullopt, std::nullopt, true);
}

void TPEvalEngineViaParsimony::Initialize() {
  // Set all PVs to Zero
  ZeroPVs();
  // Populate Leaves with Site Patterns.
  PopulateLeafParsimonyPVsWithSitePatterns();
  // Populate rootward and leafward PVs.
  PopulatePVs();
  // Compute Scores.
  ComputeScores();
}

void TPEvalEngineViaParsimony::ZeroPVs() {
  for (EdgeId edge_id = 0; edge_id < GetPVs().GetCount(); edge_id++) {
    for (const auto pv_type : PSVTypeEnum::Iterator()) {
      GetPVs().GetPV(pv_type, edge_id).setZero();
    }
  }
}

void TPEvalEngineViaParsimony::PopulatePVs() {
  // Rootward Pass (populate P PVs)
  PopulateRootwardPVs();
  // Leafward Pass (populate R PVs)
  PopulateLeafwardPVs();
}

void TPEvalEngineViaParsimony::PopulateRootwardPVs() {
  const auto rootward_node_ids = GetDAG().RootwardNodeTraversalTrace(false);
  for (const auto node_id : rootward_node_ids) {
    PopulateRootwardParsimonyPVForNode(node_id);
  }
}

void TPEvalEngineViaParsimony::PopulateLeafwardPVs() {
  const auto leafward_node_ids = GetDAG().LeafwardNodeTraversalTrace(true);
  for (const auto node_id : leafward_node_ids) {
    PopulateLeafwardParsimonyPVForNode(node_id);
  }
}

void TPEvalEngineViaParsimony::GrowNodeData(
    const size_t node_count, std::optional<const Reindexer> node_reindexer,
    std::optional<const size_t> explicit_alloc, const bool on_init) {
  // No node data to resize.
}

void TPEvalEngineViaParsimony::GrowEdgeData(
    const size_t edge_count, std::optional<const Reindexer> edge_reindexer,
    std::optional<const size_t> explicit_alloc, const bool on_init) {
  // Build resizer for resizing data.
  Resizer resizer =
      Resizer(GetTPEngine().GetEdgeCount(), GetTPEngine().GetSpareEdgeCount(),
              GetTPEngine().GetAllocatedEdgeCount(), edge_count, std::nullopt,
              explicit_alloc, GetTPEngine().GetResizingFactor());
  resizer.ApplyResizeToEigenVector<EigenVectorXd, double>(GetTopTreeScores(),
                                                          DOUBLE_NEG_INF);

  GetPVs().Resize(resizer.GetNewCount(), resizer.GetNewAlloc(), resizer.GetNewSpare());
  // Reindex work space to realign with DAG.
  if (edge_reindexer.has_value()) {
    auto pv_reindexer = GetPVs().BuildPVReindexer(
        edge_reindexer.value(), resizer.GetOldCount(), resizer.GetNewCount());
    GetPVs().Reindex(pv_reindexer);
    Reindexer::ReindexInPlace<EigenVectorXd, double>(
        GetTopTreeScores(), edge_reindexer.value(), resizer.GetNewCount());
  }
}

void TPEvalEngineViaParsimony::GrowSpareNodeData(const size_t new_node_spare_count) {
  if (new_node_spare_count > GetTPEngine().GetSpareNodeCount()) {
    GetTPEngine().GrowSpareNodeData(new_node_spare_count);
  }
  GrowNodeData(GetTPEngine().GetNodeCount());
}

void TPEvalEngineViaParsimony::GrowSpareEdgeData(const size_t new_edge_spare_count) {
  if (new_edge_spare_count > GetTPEngine().GetSpareEdgeCount()) {
    GetTPEngine().GrowSpareEdgeData(new_edge_spare_count);
  }
  GrowEdgeData(GetTPEngine().GetEdgeCount());
}

void TPEvalEngineViaParsimony::GrowEngineForDAG(
    std::optional<Reindexer> node_reindexer, std::optional<Reindexer> edge_reindexer) {
  const auto &dag = GetTPEngine().GetDAG();
  GrowNodeData(dag.NodeCount(), node_reindexer, std::nullopt, false);
  GrowEdgeData(dag.EdgeCountWithLeafSubsplits(), edge_reindexer, std::nullopt, false);
}

void TPEvalEngineViaParsimony::GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                                         const bool via_reference,
                                                         const bool use_unique_temps) {
  if (via_reference) {
    if (use_unique_temps) {
      GrowSpareNodeData(spare_nodes_per_nni_ * adjacent_nnis.size());
      GrowSpareEdgeData(spare_edges_per_nni_ * adjacent_nnis.size());
    } else {
      GrowSpareNodeData(spare_nodes_per_nni_);
      GrowSpareEdgeData(spare_edges_per_nni_);
    }
    GrowSpareEdgeData(adjacent_nnis.size());
  } else {
    GrowNodeData(GetGraftDAG().NodeCountWithoutDAGRoot());
    GrowEdgeData(GetGraftDAG().EdgeCountWithLeafSubsplits());
  }
}

void TPEvalEngineViaParsimony::UpdateEngineAfterDAGAddNodePair(
    const NNIOperation &post_nni, const NNIOperation &pre_nni,
    std::optional<size_t> new_tree_id) {
  // Copy over edge data.
  GetTPEngine().CopyOverEdgeDataFromPreNNIToPostNNI(
      post_nni, pre_nni,
      [this](const EdgeId src, const EdgeId dest) { CopyEdgeData(src, dest); },
      new_tree_id);
  // Populate PVs.
  PopulatePVs();
}

void TPEvalEngineViaParsimony::UpdateEngineAfterModifyingDAG(
    const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
    const size_t prev_node_count, const Reindexer &node_reindexer,
    const size_t prev_edge_count, const Reindexer &edge_reindexer) {
  Initialize();
  // Update scores.
  ComputeScores();
}

double TPEvalEngineViaParsimony::GetTopTreeScoreWithEdge(const EdgeId edge_id) const {
  return TPEvalEngine::GetTopTreeScoreWithEdge(edge_id);
}

double TPEvalEngineViaParsimony::GetTopTreeScoreWithProposedNNI(
    const NNIOperation &post_nni, const NNIOperation &pre_nni,
    const size_t spare_offset, std::optional<BitsetEdgeIdMap> best_edge_map) {
  using NNIClade = NNIOperation::NNIClade;
  using NNICladeEnum = NNIOperation::NNICladeEnum;
  GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits());
  // Node ids from pre-NNI in DAG.
  NNICladeEnum::Array<EdgeId> pre_id_map;
  NNICladeEnum::Array<EdgeId> post_id_map;
  const auto pre_edge_id = GetDAG().GetEdgeIdx(pre_nni);
  // Create mapping between pre-NNI and post-NNI.
  const auto clade_map =
      NNIOperation::BuildNNICladeMapFromPreNNIToNNI(pre_nni, post_nni);
  // PLV ids from pre-NNI in DAG.
  auto choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(pre_edge_id);
  pre_id_map[NNIClade::ParentFocal] = choices.parent;
  pre_id_map[NNIClade::ParentSister] = choices.sister;
  pre_id_map[NNIClade::ChildLeft] = choices.left_child;
  pre_id_map[NNIClade::ChildRight] = choices.right_child;
  // Use clade mapping to reference pre-NNI PVs for post-NNI PVs.
  for (const auto nni_clade : NNICladeEnum::Iterator()) {
    post_id_map[nni_clade] = pre_id_map[clade_map[nni_clade]];
  }
  // Get temp PVs for post-NNI PVs and edge lengths.
  const PVId q_pvid = PVId(spare_offset * 3);
  const PVId pleft_pvid = PVId((spare_offset * 3) + 1);
  const PVId pright_pvid = PVId((spare_offset * 3) + 2);
  // Compute Pleft and Pright.
  for (size_t pattern_idx = 0; pattern_idx < GetSitePattern().PatternCount();
       pattern_idx++) {
    GetPVs().GetPV(pleft_pvid).col(pattern_idx) =
        ParentPartial(TotalPPartial(post_id_map[NNIClade::ChildLeft], pattern_idx));
    GetPVs().GetPV(pright_pvid).col(pattern_idx) =
        ParentPartial(TotalPPartial(post_id_map[NNIClade::ChildRight], pattern_idx));
  }
  // Compute Q.
  for (size_t pattern_idx = 0; pattern_idx < GetSitePattern().PatternCount();
       pattern_idx++) {
    auto partials_from_parent =
        ParentPartial(GetPVs()
                          .GetPV(PSVType::Q, post_id_map[NNIClade::ParentFocal])
                          .col(pattern_idx));
    for (const auto child_id :
         {post_id_map[NNIClade::ChildLeft], post_id_map[NNIClade::ChildRight]}) {
      EdgeId sister_id = ((child_id == post_id_map[NNIClade::ChildLeft])
                              ? post_id_map[NNIClade::ChildRight]
                              : post_id_map[NNIClade::ChildLeft]);
      auto partials_from_sister = ParentPartial(TotalPPartial(sister_id, pattern_idx));
      GetPVs().GetPV(q_pvid).col(pattern_idx) =
          partials_from_sister + partials_from_parent;
    }
  }
  // Compute total parsimony.
  double score = ParsimonyScore(q_pvid, pleft_pvid, pright_pvid);
  return score;
}

void TPEvalEngineViaParsimony::CopyEdgeData(const EdgeId src_edge_id,
                                            const EdgeId dest_edge_id) {
  TPEvalEngine::CopyEdgeData(src_edge_id, dest_edge_id);
}

void TPEvalEngineViaParsimony::PopulateLeafParsimonyPVsWithSitePatterns() {
  // first check that the psv_handler has been resized to deal with the leaf labels
  Assert(GetPVs().GetCount() >= GetSitePattern().TaxonCount(),
         "Error in SankoffHandler::GenerateLeafPartials: "
         "parsimony_pvs_ should be initialized to accomodate"
         "the number of leaf nodes in the GetSitePattern().");

  // Iterate over all leaf nodes to instantiate each with P partial values
  for (const auto node_id : GetDAG().GetLeafNodeIds()) {
    for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
      for (const auto adj_node_id :
           GetDAG().GetDAGNode(node_id).GetNeighbors(Direction::Rootward, clade)) {
        const auto edge_id = GetDAG().GetEdgeIdx(adj_node_id, node_id);
        SankoffPartial leaf_partials(state_count_, GetSitePattern().PatternCount());
        // set leaf node partial to have big_double_ infinity substitute
        leaf_partials.block(0, 0, state_count_, GetSitePattern().PatternCount())
            .fill(big_double_);
        // now fill in appropriate entries of the leaf-partial where non-infinite
        for (size_t pattern_idx = 0; pattern_idx < GetSitePattern().PatternCount();
             pattern_idx++) {
          auto site_val =
              GetSitePattern().GetPatternSymbol(node_id.value_, pattern_idx);
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
        GetPVs().GetPV(PSVType::PLeft, edge_id) = leaf_partials;
        GetPVs().GetPV(PSVType::PRight, edge_id).fill(0);
      }
    }
  }
}

void TPEvalEngineViaParsimony::PopulateRootwardParsimonyPVForNode(
    const NodeId node_id) {
  for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto adj_node_id :
         GetDAG().GetDAGNode(node_id).GetNeighbors(Direction::Rootward, clade)) {
      const EdgeId edge_id = GetDAG().GetEdgeIdx(adj_node_id, node_id);
      // Populate edge PLV by accumulating parsimony up given edge.
      PopulateRootwardParsimonyPVForEdge(edge_id);
    }
  }
}

void TPEvalEngineViaParsimony::PopulateRootwardParsimonyPVForEdge(
    const EdgeId edge_id) {
  const auto choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(edge_id);
  if (choices.left_child != NoId && choices.right_child != NoId) {
    // Accumulate parsimony from left and right child.
    const EdgeId left_child_edge_id = choices.left_child;
    const EdgeId right_child_edge_id = choices.right_child;
    PopulateRootwardParsimonyPVForEdge(edge_id, left_child_edge_id,
                                       right_child_edge_id);
  }
}

void TPEvalEngineViaParsimony::PopulateLeafwardParsimonyPVForNode(
    const NodeId node_id) {
  for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto adj_node_id :
         GetDAG().GetDAGNode(node_id).GetNeighbors(Direction::Leafward, clade)) {
      const EdgeId edge_id = GetDAG().GetEdgeIdx(node_id, adj_node_id);
      // Populate edge PLV by accumulating parsimony down given edge.
      PopulateLeafwardParsimonyPVForEdge(edge_id);
    }
  }
}

void TPEvalEngineViaParsimony::PopulateLeafwardParsimonyPVForEdge(
    const EdgeId edge_id) {
  const auto choices = GetTPEngine().GetChoiceMap().GetEdgeChoice(edge_id);
  if (choices.parent != NoId && choices.sister != NoId) {
    // Evolve down parent.
    const auto edge = GetDAG().GetDAGEdge(edge_id);
    const EdgeId parent_edge_id = choices.parent;
    const EdgeId sister_edge_id = choices.sister;
    const EdgeId left_child_edge_id =
        (edge.GetSubsplitClade() == SubsplitClade::Left) ? edge_id : sister_edge_id;
    const EdgeId right_child_edge_id =
        (edge.GetSubsplitClade() == SubsplitClade::Left) ? sister_edge_id : edge_id;

    PopulateLeafwardParsimonyPVForEdge(parent_edge_id, left_child_edge_id,
                                       right_child_edge_id);
  }
}

void TPEvalEngineViaParsimony::ComputeScores(std::optional<EdgeIdVector> opt_edge_ids) {
  for (EdgeId edge_id = 0; edge_id < GetDAG().EdgeCountWithLeafSubsplits(); edge_id++) {
    GetTopTreeScores()[edge_id.value_] = ParsimonyScore(edge_id);
  }
}

EigenVectorXd TPEvalEngineViaParsimony::ParentPartial(EigenVectorXd child_partials) {
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

EigenVectorXd TPEvalEngineViaParsimony::TotalPPartial(const EdgeId edge_id,
                                                      const size_t site_idx) {
  return TotalPPartial(GetPVs().GetPVIndex(PSVType::PLeft, edge_id),
                       GetPVs().GetPVIndex(PSVType::PRight, edge_id), site_idx);
}

EigenVectorXd TPEvalEngineViaParsimony::TotalPPartial(const PVId edge_pleft_pvid,
                                                      const PVId edge_pright_pvid,
                                                      const size_t site_idx) {
  return GetPVs().GetPV(edge_pleft_pvid).col(site_idx) +
         GetPVs().GetPV(edge_pright_pvid).col(site_idx);
}

void TPEvalEngineViaParsimony::PopulateRootwardParsimonyPVForEdge(
    const EdgeId parent_id, const EdgeId left_child_id, const EdgeId right_child_id) {
  for (size_t pattern_idx = 0; pattern_idx < GetSitePattern().PatternCount();
       pattern_idx++) {
    // Which child partial is in right or left doesn't actually matter because they
    // are summed when calculating q_partials.
    GetPVs().GetPV(PSVType::PLeft, parent_id).col(pattern_idx) =
        ParentPartial(TotalPPartial(left_child_id, pattern_idx));
    GetPVs().GetPV(PSVType::PRight, parent_id).col(pattern_idx) =
        ParentPartial(TotalPPartial(right_child_id, pattern_idx));
  }
}

void TPEvalEngineViaParsimony::PopulateLeafwardParsimonyPVForEdge(
    const EdgeId parent_id, const EdgeId left_child_id, const EdgeId right_child_id) {
  for (size_t pattern_idx = 0; pattern_idx < GetSitePattern().PatternCount();
       pattern_idx++) {
    auto partials_from_parent =
        ParentPartial(GetPVs().GetPV(PSVType::Q, parent_id).col(pattern_idx));
    for (const auto child_id : {left_child_id, right_child_id}) {
      EdgeId sister_id = ((child_id == left_child_id) ? right_child_id : left_child_id);
      auto partials_from_sister = ParentPartial(TotalPPartial(sister_id, pattern_idx));
      GetPVs().GetPV(PSVType::Q, child_id).col(pattern_idx) =
          partials_from_sister + partials_from_parent;
    }
  }
}

double TPEvalEngineViaParsimony::ParsimonyScore(const EdgeId edge_id) {
  return ParsimonyScore(GetPVs().GetPVIndex(PSVType::Q, edge_id),
                        GetPVs().GetPVIndex(PSVType::PLeft, edge_id),
                        GetPVs().GetPVIndex(PSVType::PRight, edge_id));
}

double TPEvalEngineViaParsimony::ParsimonyScore(const PVId edge_q_pvid,
                                                const PVId edge_pleft_pvid,
                                                const PVId edge_pright_pvid) {
  auto weights = GetSitePattern().GetWeights();
  double total_parsimony = 0.;
  for (size_t pattern = 0; pattern < GetSitePattern().PatternCount(); pattern++) {
    // Note: doing ParentPartial first for the left and right p_partials and then
    // adding them together will give the same minimum parsimony score, but doesn't
    // give correct Sankoff Partial vector for the new rooting
    auto total_tree =
        ParentPartial(TotalPPartial(edge_pleft_pvid, edge_pright_pvid, pattern));
    total_tree += ParentPartial(GetPVs().GetPV(edge_q_pvid).col(pattern));

    // If node_id is the root node, calculating the total_tree vector like so does not
    // yield the SankoffPartial of an actual rooting, but this will not change the
    // minimum value in the partial, so the root node can still be used to calculate
    // the parsimony score.
    total_parsimony +=
        *std::min_element(total_tree.begin(), total_tree.end()) * weights[pattern];
  }
  return total_parsimony;
}
