// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "nni_engine.hpp"
#include "stopwatch.hpp"

using PLVType = PLVNodeHandler::PLVType;

NNIEngine::NNIEngine(GPDAG &dag, std::optional<GPEngine *> gp_engine,
                     std::optional<TPEngine *> tp_engine)
    : dag_(dag), graft_dag_(std::make_unique<GraftDAG>(dag)) {
  if (gp_engine.has_value() && gp_engine.value()) {
    MakeGPEvalEngine(gp_engine.value());
  }
  if (tp_engine.has_value() && tp_engine.value()) {
    MakeTPEvalEngine(tp_engine.value());
  }
}

// ** Access

double NNIEngine::GetScoreByNNI(const NNIOperation &nni) const {
  // Check if score has already been stored in NNI Engine.
  if (GetScoredNNIs().find(nni) != GetScoredNNIs().end()) {
    return GetScoredNNIs().at(nni);
  }
  if (GetPastScoredNNIs().find(nni) != GetScoredNNIs().end()) {
    return GetPastScoredNNIs().at(nni);
  }
  // Otherwise, check if score is in NNI Eval Engine.
  return GetEvalEngine().ScoreInternalNNIByNNI(nni);
}

double NNIEngine::GetScoreByEdge(const EdgeId edge_id) const {
  // Check for score in NNI Eval Engine.
  return GetEvalEngine().ScoreInternalNNIByEdge(edge_id);
}

// ** NNI Evaluation Engine

NNIEvalEngineViaGP &NNIEngine::MakeGPEvalEngine(GPEngine *gp_engine) {
  Assert(gp_engine != nullptr, "Cannot MakeGPEvalEngine with nullptr.");
  eval_engine_via_gp_ = std::make_unique<NNIEvalEngineViaGP>(*this, *gp_engine);
  SelectGPEvalEngine();
  return GetGPEvalEngine();
}

NNIEvalEngineViaTP &NNIEngine::MakeTPEvalEngine(TPEngine *tp_engine) {
  Assert(tp_engine != nullptr, "Cannot MakeTPEvalEngine with nullptr.");
  eval_engine_via_tp_ = std::make_unique<NNIEvalEngineViaTP>(*this, *tp_engine);
  SelectTPLikelihoodEvalEngine();
  return GetTPEvalEngine();
}

void NNIEngine::ClearEvalEngineInUse() {
  for (auto eval_engine_type : NNIEvalEngineTypeEnum::Iterator()) {
    eval_engine_in_use_[eval_engine_type] = false;
  }
}

void NNIEngine::SelectEvalEngine(const NNIEvalEngineType eval_engine_type) {
  switch (eval_engine_type) {
    case NNIEvalEngineType::GPEvalEngine:
      SelectGPEvalEngine();
      break;
    case NNIEvalEngineType::TPEvalEngineViaLikelihood:
      SelectTPLikelihoodEvalEngine();
      break;
    case NNIEvalEngineType::TPEvalEngineViaParsimony:
      SelectTPParsimonyEvalEngine();
      break;
    default:
      Failwith("Invalid NNIEvalEngineType.");
  }
}

void NNIEngine::SelectGPEvalEngine() {
  Assert(HasGPEvalEngine(), "Must MakeGPEvalEngine before selecting it.");
  ClearEvalEngineInUse();
  eval_engine_in_use_[NNIEvalEngineType::GPEvalEngine] = true;
  eval_engine_ = &GetGPEvalEngine();
}

void NNIEngine::SelectTPLikelihoodEvalEngine() {
  Assert(HasTPEvalEngine() && GetTPEngine().HasLikelihoodEvalEngine(),
         "Must MakeTPEvalEngine with LikelihoodEvalEngine before selecting it.");
  ClearEvalEngineInUse();
  eval_engine_in_use_[NNIEvalEngineType::TPEvalEngineViaLikelihood] = true;
  GetTPEngine().SelectLikelihoodEvalEngine();
  eval_engine_ = &GetTPEvalEngine();
}

void NNIEngine::SelectTPParsimonyEvalEngine() {
  Assert(HasTPEvalEngine() && GetTPEngine().HasParsimonyEvalEngine(),
         "Must MakeTPEvalEngine with ParsimonyEvalEngine before selecting it.");
  ClearEvalEngineInUse();
  eval_engine_in_use_[NNIEvalEngineType::TPEvalEngineViaParsimony];
  GetTPEngine().SelectParsimonyEvalEngine();
  eval_engine_ = &GetTPEvalEngine();
}

void NNIEngine::InitEvalEngine() {
  if (IsEvalEngineInUse(NNIEvalEngineType::GPEvalEngine)) {
    GetGPEvalEngine().Init();
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaLikelihood)) {
    GetTPEvalEngine().Init();
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaParsimony)) {
    GetTPEvalEngine().Init();
  }
}

void NNIEngine::PrepEvalEngine() {
  if (IsEvalEngineInUse(NNIEvalEngineType::GPEvalEngine)) {
    GetGPEvalEngine().Prep();
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaLikelihood)) {
    GetTPEvalEngine().Prep();
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaParsimony)) {
    GetTPEvalEngine().Prep();
  }
}

void NNIEngine::GrowEvalEngineForDAG(std::optional<Reindexer> node_reindexer,
                                     std::optional<Reindexer> edge_reindexer) {
  if (IsEvalEngineInUse(NNIEvalEngineType::GPEvalEngine)) {
    GetGPEvalEngine().GrowEngineForDAG(node_reindexer, edge_reindexer);
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaLikelihood)) {
    GetTPEvalEngine().GrowEngineForDAG(node_reindexer, edge_reindexer);
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaParsimony)) {
    GetTPEvalEngine().GrowEngineForDAG(node_reindexer, edge_reindexer);
  }
}

void NNIEngine::GrowEvalEngineForAdjacentNNIs(const bool via_reference,
                                              const bool use_unique_temps) {
  if (IsEvalEngineInUse(NNIEvalEngineType::GPEvalEngine)) {
    GetGPEvalEngine().GrowEngineForAdjacentNNIs(GetAdjacentNNIs(), via_reference,
                                                use_unique_temps);
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaLikelihood)) {
    GetTPEvalEngine().GrowEngineForAdjacentNNIs(GetAdjacentNNIs(), via_reference,
                                                use_unique_temps);
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaParsimony)) {
    GetTPEvalEngine().GrowEngineForAdjacentNNIs(GetAdjacentNNIs(), via_reference,
                                                use_unique_temps);
  }
}

void NNIEngine::UpdateEvalEngineAfterModifyingDAG(
    const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
    const size_t prev_node_count, const Reindexer &node_reindexer,
    const size_t prev_edge_count, const Reindexer &edge_reindexer) {
  if (IsEvalEngineInUse(NNIEvalEngineType::GPEvalEngine)) {
    auto node_reindexer_without_dag(node_reindexer);
    node_reindexer_without_dag =
        node_reindexer.RemoveNewIndex(GetDAG().GetDAGRootNodeId().value_);
    GetGPEvalEngine().UpdateEngineAfterModifyingDAG(nni_to_pre_nni, prev_node_count,
                                                    node_reindexer_without_dag,
                                                    prev_edge_count, edge_reindexer);
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaLikelihood)) {
    GetTPEvalEngine().UpdateEngineAfterModifyingDAG(nni_to_pre_nni, prev_node_count,
                                                    node_reindexer, prev_edge_count,
                                                    edge_reindexer);
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaParsimony)) {
    GetTPEvalEngine().UpdateEngineAfterModifyingDAG(nni_to_pre_nni, prev_node_count,
                                                    node_reindexer, prev_edge_count,
                                                    edge_reindexer);
  }
}

void NNIEngine::ScoreAdjacentNNIs() {
  if (IsEvalEngineInUse(NNIEvalEngineType::GPEvalEngine)) {
    GetGPEvalEngine().ScoreAdjacentNNIs(GetAdjacentNNIs());
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaLikelihood)) {
    GetTPEvalEngine().ScoreAdjacentNNIs(GetAdjacentNNIs());
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaParsimony)) {
    GetTPEvalEngine().ScoreAdjacentNNIs(GetAdjacentNNIs());
  }
}

// ** Runners

void NNIEngine::Run(const bool is_quiet) {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cout);
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);
  os << "RunInit(): ";
  // Initialize Once at the start of the loop
  RunInit(is_quiet);
  os << timer.Lap() << std::endl;
  // Loop until no more eligible NNIs are found.
  while (GetAdjacentNNICount() > 0) {
    os << "RunMainLoop(): ";
    RunMainLoop(is_quiet);
    os << timer.Lap() << std::endl;

    os << "RunPostLoop(): ";
    RunPostLoop(is_quiet);
    os << timer.Lap() << std::endl;
  }
}

void NNIEngine::RunInit(const bool is_quiet) {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cout);
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);
  // Initialize Adjacent NNIs based on starting state of DAG.
  ResetAllNNIs();
  os << "RunInit::ResetAllNNIs: " << timer.Lap() << std::endl;
  SyncAdjacentNNIsWithDAG();
  os << "RunInit::SyncAdjacentNNIsWithDAG: " << timer.Lap() << std::endl;
  PrepEvalEngine();
  os << "RunInit::PrepEvalEngine: " << timer.Lap() << std::endl;
  FilterInit();
  os << "RunInit::FilterInit: " << timer.Lap() << std::endl;
}

void NNIEngine::RunMainLoop(const bool is_quiet) {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cout);
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);

  // (1) Add all adjacent NNIs to the GraftDAG.
  GraftAdjacentNNIsToDAG();
  os << "RunMainLoop::GraftAdjacentNNIsToDAG: " << timer.Lap() << std::endl;
  // (2) Evaluate each adjacent NNI.
  FilterPreUpdate();
  os << "RunMainLoop::FilterPreUpdate: " << timer.Lap() << std::endl;
  FilterEvaluateAdjacentNNIs();
  os << "RunMainLoop::FilterEvaluateAdjacentNNIs: " << timer.Lap() << std::endl;
  FilterPostUpdate();
  os << "RunMainLoop::FilterPostUpdate: " << timer.Lap() << std::endl;
  // (3) Select whether to accept or reject adjacent NNIs via filter.
  FilterProcessAdjacentNNIs();
  os << "RunMainLoop::FilterProcess: " << timer.Lap() << std::endl;
  // (4a) Remove adjacent NNIs from GraftDAG.
  RemoveAllGraftedNNIsFromDAG();
  // (4b) Add accepted NNIs permanently to DAG.
  AddAcceptedNNIsToDAG(is_quiet);
  os << "RunMainLoop::RemoveAndAddNNIs: " << timer.Lap() << std::endl;

  iter_count_++;
}

void NNIEngine::RunPostLoop(const bool is_quiet) {
  // (5a) Update Adjacent NNIs to reflect added NNI.
  UpdateAdjacentNNIs(true);
  // (5b) Reset Accepted NNIs and save results.
  UpdateAcceptedNNIs(true);
  // (5c) Reset Rejected NNIs and save results.
  UpdateRejectedNNIs(true);
  // (5d) Reset Scored NNIs and save results.
  UpdateScoredNNIs(true);
}

// ** Filter Functions

void NNIEngine::SetNoEvaluate() {
  SetFilterEvalFunction([](NNIEngine &this_nni_engine, NNIEvalEngine &this_eval_engine,
                           GraftDAG &this_graft_dag,
                           const NNIOperation &nni) -> double { return 0.0; });
}

void NNIEngine::SetNoFilter(const bool set_nni_to_pass) {
  SetFilterProcessFunction(
      [set_nni_to_pass](NNIEngine &this_nni_engine, NNIEvalEngine &this_eval_engine,
                        GraftDAG &this_graft_dag, const NNIOperation &nni,
                        const double nni_score) -> bool { return set_nni_to_pass; });
}

void NNIEngine::SetFilterBySetOfNNIs(const std::set<NNIOperation> &nnis_to_accept) {
  SetFilterProcessFunction(
      [&nnis_to_accept](NNIEngine &this_nni_engine, NNIEvalEngine &this_eval_engine,
                        GraftDAG &this_graft_dag, const NNIOperation &nni,
                        const double nni_score) -> bool {
        bool nni_passes = (nnis_to_accept.find(nni) != nnis_to_accept.end());
        return nni_passes;
      });
}

void NNIEngine::SetMinScoreCutoff(const double score_cutoff) {
  SetFilterProcessFunction(
      [score_cutoff](NNIEngine &this_nni_engine, NNIEvalEngine &this_eval_engine,
                     GraftDAG &this_graft_dag, const NNIOperation &nni,
                     const double nni_score) -> bool {
        bool nni_passes = (nni_score >= score_cutoff);
        return nni_passes;
      });
}

void NNIEngine::SetMaxScoreCutoff(const double score_cutoff) {
  SetFilterProcessFunction(
      [score_cutoff](NNIEngine &this_nni_engine, NNIEvalEngine &this_eval_engine,
                     GraftDAG &this_graft_dag, const NNIOperation &nni,
                     const double nni_score) -> bool {
        bool nni_passes = (nni_score <= score_cutoff);
        return nni_passes;
      });
}

void NNIEngine::SetScoredNNIsFromEvalEngine() {
  // SetFilterPostUpdateFunction([](NNIEngine &this_nni_engine,
  //                                NNIEvalEngine &this_eval_engine,
  //                                GraftDAG &this_graft_dag) {
  //   this_nni_engine.GetScoredNNIs().clear();
  //   for (const auto &[nni, score] : this_eval_engine.GetScoredNNIs()) {
  //     if (!this_graft_dag.ContainsNNI(nni)) {
  //       this_nni_engine.GetScoredNNIs()[nni] = score;
  //     }
  //   }
  // });
}

// ** Filter Subroutines

void NNIEngine::FilterInit() {
  if (filter_init_fn_) {
    (filter_init_fn_)(*this, GetEvalEngine(), GetGraftDAG());
  }
}

void NNIEngine::FilterPreUpdate() {
  if (filter_pre_update_fn_) {
    (filter_pre_update_fn_)(*this, GetEvalEngine(), GetGraftDAG());
  }
}

void NNIEngine::FilterEvaluateAdjacentNNIs() {
  if (!filter_eval_fn_) return;
  for (const auto &nni : GetAdjacentNNIs()) {
    const double nni_score =
        (filter_eval_fn_)(*this, GetEvalEngine(), GetGraftDAG(), nni);
    AddScoreForNNI(nni, nni_score);
  }
}

void NNIEngine::FilterPostUpdate() {
  if (filter_post_update_fn_) {
    (filter_post_update_fn_)(*this, GetEvalEngine(), GetGraftDAG());
  }
}

void NNIEngine::FilterProcessAdjacentNNIs() {
  Assert(filter_process_fn_, "Must assign a filter process function.");
  for (const auto &nni : GetAdjacentNNIs()) {
    double nni_score = (*GetScoredNNIs().find(nni)).second;
    const bool accept_nni =
        (filter_process_fn_)(*this, GetEvalEngine(), GetGraftDAG(), nni, nni_score);
    if (accept_nni) {
      accepted_nnis_.insert(nni);
    } else {
      rejected_nnis_.insert(nni);
    }
  }
}

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

// ** Filtering Scheme

void NNIEngine::SetGPLikelihoodCutoffFilteringScheme(const double score_cutoff) {
  Assert(HasGPEvalEngine(), "Must MakeGPEvalEngine before using the filtering scheme.");
  SelectGPEvalEngine();
  SetFilterPreUpdateFunction(
      [](NNIEngine &this_nni_engine, NNIEvalEngine &this_eval_engine,
         GraftDAG &this_graft_dag) { this_nni_engine.ScoreAdjacentNNIs(); });
  SetScoredNNIsFromEvalEngine();
  SetMinScoreCutoff(score_cutoff);
}

void NNIEngine::SetTPLikelihoodCutoffFilteringScheme(const double score_cutoff) {
  Assert(HasTPEvalEngine() && GetTPEngine().HasLikelihoodEvalEngine(),
         "Must MakeTPEvalEngine before using the filtering scheme.");
  SelectTPLikelihoodEvalEngine();
  SetFilterPreUpdateFunction(
      [](NNIEngine &this_nni_engine, NNIEvalEngine &this_eval_engine,
         GraftDAG &this_graft_dag) { this_nni_engine.ScoreAdjacentNNIs(); });
  SetScoredNNIsFromEvalEngine();
  SetMinScoreCutoff(score_cutoff);
}

void NNIEngine::SetTPParsimonyCutoffFilteringScheme(const double score_cutoff) {
  Assert(HasTPEvalEngine() && GetTPEngine().HasParsimonyEvalEngine(),
         "Must MakeTPEvalEngine before using the filtering scheme.");
  SelectTPParsimonyEvalEngine();
  SetFilterPreUpdateFunction(
      [](NNIEngine &this_nni_engine, NNIEvalEngine &this_eval_engine,
         GraftDAG &this_graft_dag) { this_nni_engine.ScoreAdjacentNNIs(); });
  SetScoredNNIsFromEvalEngine();
  SetMaxScoreCutoff(score_cutoff);
}

void NNIEngine::SetGPLikelihoodDropFilteringScheme(const double score_cutoff) {
  Assert(HasGPEvalEngine(), "Must MakeGPEvalEngine before using the filtering scheme.");
  SelectGPEvalEngine();
  SetFilterPreUpdateFunction([score_cutoff](NNIEngine &this_nni_engine,
                                            NNIEvalEngine &this_eval_engine,
                                            GraftDAG &this_graft_dag) {
    this_nni_engine.ScoreAdjacentNNIs();
    double max = this_eval_engine.GetMaxScore();
    this_nni_engine.SetMinScoreCutoff(max - score_cutoff);
  });
  SetScoredNNIsFromEvalEngine();
}

void NNIEngine::SetTPLikelihoodDropFilteringScheme(const double score_cutoff) {
  Assert(HasTPEvalEngine() && GetTPEngine().HasLikelihoodEvalEngine(),
         "Must MakeTPEvalEngine before using the filtering scheme.");
  SelectTPLikelihoodEvalEngine();
  SetFilterPreUpdateFunction([score_cutoff](NNIEngine &this_nni_engine,
                                            NNIEvalEngine &this_eval_engine,
                                            GraftDAG &this_graft_dag) {
    this_nni_engine.ScoreAdjacentNNIs();
    double max = this_eval_engine.GetMaxScore();
    this_nni_engine.SetMinScoreCutoff(max - score_cutoff);
  });
  SetScoredNNIsFromEvalEngine();
}

void NNIEngine::SetTPParsimonyDropFilteringScheme(const double score_cutoff) {
  Assert(HasTPEvalEngine() && GetTPEngine().HasParsimonyEvalEngine(),
         "Must MakeTPEvalEngine before using the filtering scheme.");
  SelectTPParsimonyEvalEngine();
  SetFilterPreUpdateFunction([score_cutoff](NNIEngine &this_nni_engine,
                                            NNIEvalEngine &this_eval_engine,
                                            GraftDAG &this_graft_dag) {
    this_nni_engine.ScoreAdjacentNNIs();
    double min = this_eval_engine.GetMinScore();
    this_nni_engine.SetMaxScoreCutoff(min + score_cutoff);
  });
  SetScoredNNIsFromEvalEngine();
}

void NNIEngine::SetTopNScoreFilteringScheme(const size_t n, const bool max_is_best) {
  SetFilterPreUpdateFunction([n, max_is_best](NNIEngine &this_nni_engine,
                                              NNIEvalEngine &this_eval_engine,
                                              GraftDAG &this_graft_dag) {
    this_nni_engine.ScoreAdjacentNNIs();
    DoubleVector nni_scores;
    for (const auto &[nni, score] : this_eval_engine.GetScoredNNIs()) {
      std::ignore = nni;
      nni_scores.push_back(score);
    }
    double score_cutoff = -INFINITY;
    if (!nni_scores.empty()) {
      std::sort(nni_scores.begin(), nni_scores.end());
      if (max_is_best) {
        score_cutoff = nni_scores[nni_scores.size() - std::min(n, nni_scores.size())];
      } else {
        score_cutoff = nni_scores[std::min(n, nni_scores.size()) - 1];
      }
    }
    if (max_is_best) {
      this_nni_engine.SetMinScoreCutoff(score_cutoff);
    } else {
      this_nni_engine.SetMaxScoreCutoff(score_cutoff);
    }
    this_nni_engine.SetScoredNNIsFromEvalEngine();
  });
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
    const NNIEngine::KeyIndexMap &pre_key_idx) const {
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
      PLVNodeHandler::GetPVIndex(PLVType::RHat, parent_id, node_count).value_;
  key_idx_map[KeyIndex::Parent_RFocal] =
      PLVNodeHandler::GetPVIndex(PLVTypeEnum::RPLVType(!is_left_clade_sister),
                                 parent_id, node_count)
          .value_;
  key_idx_map[KeyIndex::Parent_PHatSister] =
      PLVNodeHandler::GetPVIndex(PLVTypeEnum::PPLVType(is_left_clade_sister), parent_id,
                                 node_count)
          .value_;
  key_idx_map[KeyIndex::Child_P] =
      PLVNodeHandler::GetPVIndex(PLVType::P, child_id, node_count).value_;
  key_idx_map[KeyIndex::Child_PHatLeft] =
      PLVNodeHandler::GetPVIndex(PLVType::PHatLeft, child_id, node_count).value_;
  key_idx_map[KeyIndex::Child_PHatRight] =
      PLVNodeHandler::GetPVIndex(PLVType::PHatRight, child_id, node_count).value_;

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

// ** DAG Maintenance

void NNIEngine::AddAcceptedNNIsToDAG(const bool is_quiet) {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cout);
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);
  // Initialize reindexers for remapping after adding nodes.
  const size_t prev_node_count = GetDAG().NodeCount();
  const size_t prev_edge_count = GetDAG().EdgeCountWithLeafSubsplits();
  node_reindexer_ = Reindexer::IdentityReindexer(GetDAG().NodeCount());
  edge_reindexer_ = Reindexer::IdentityReindexer(GetDAG().EdgeCountWithLeafSubsplits());
  // Build map from NNI to pre-NNI.
  std::map<NNIOperation, NNIOperation> nni_to_pre_nni;
  for (const auto &nni : GetAcceptedNNIs()) {
    auto adj_nnis = GetDAG().FindAllNNINeighborsInDAG(nni);
    bool nni_found = false;
    for (const auto clade : SubsplitCladeEnum::Iterator()) {
      const auto adj_nni = adj_nnis[clade];
      if (adj_nni.has_value() &&
          (GetAdjacentNNIs().find(adj_nni.value()) == GetAdjacentNNIs().end())) {
        nni_to_pre_nni[nni] = adj_nni.value();
        nni_found = true;
      }
    }
    Assert(nni_found, "NNI not found to be adjacent to DAG.");
  }
  // Add NNI to DAG.
  os << "AddAcceptedNNIsToDAG: " << std::endl;
  // BitsetPairVector nodes_to_add;
  size_t nni_count = 0;
  for (const auto &nni : GetAcceptedNNIs()) {
    os << "AddAcceptedNNIsToDAG: " << nni_count++ << " of " << GetAcceptedNNIs().size()
       << std::endl;
    auto mods = GetDAG().AddNodePair(nni);
    node_reindexer_ = node_reindexer_.ComposeWith(mods.node_reindexer);
    edge_reindexer_ = edge_reindexer_.ComposeWith(mods.edge_reindexer);
    // nodes_to_add.push_back({nni.GetParent(), nni.GetChild()});
  }
  // auto mods = GetDAG().AddNodes(nodes_to_add);
  // node_reindexer_ = mods.node_reindexer;
  // edge_reindexer_ = mods.edge_reindexer;
  os << "AddAcceptedNNIsToDAG::AddAll: " << timer.Lap() << std::endl;

  GrowEvalEngineForDAG(node_reindexer_, edge_reindexer_);
  os << "AddAcceptedNNIsToDAG::GrowEvalEngine: " << timer.Lap() << std::endl;
  UpdateEvalEngineAfterModifyingDAG(nni_to_pre_nni, prev_node_count, node_reindexer_,
                                    prev_edge_count, edge_reindexer_);
  os << "AddAcceptedNNIsToDAG::UpdateEvalEngine: " << timer.Lap() << std::endl;
}

void NNIEngine::GraftAdjacentNNIsToDAG() {
  for (const auto &nni : GetAdjacentNNIs()) {
    GetGraftDAG().AddNodePair(nni);
  }
}

void NNIEngine::RemoveAllGraftedNNIsFromDAG() {
  GetGraftDAG().RemoveAllGrafts();
  graft_dag_ = std::make_unique<GraftDAG>(GetDAG());
}

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
            if (GetIncludeRootsplitNNIs() || !parent_bitset.SubsplitIsRootsplit()) {
              SafeAddOutputNNIsToAdjacentNNIs(parent_bitset, child_bitset,
                                              is_edge_on_left);
            }
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
        if (GetIncludeRootsplitNNIs() || !parent_bitset.SubsplitIsRootsplit()) {
          AddAllNNIsFromNodeVectorToAdjacentNNIs(node_id, adjacent_node_ids,
                                                 is_edge_on_left, is_edge_leafward);
        }
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
  SafeInsert(GetScoredNNIs(), nni, score);
};

void NNIEngine::UpdateAdjacentNNIs(const bool save_rejected_nnis) {
  adjacent_nnis_.clear();
  if (save_rejected_nnis) {
    adjacent_nnis_.insert(rejected_nnis_.begin(), rejected_nnis_.end());
  }
  for (const auto &nni : GetAcceptedNNIs()) {
    const auto nni_edge_id = GetDAG().GetEdgeIdx(nni);
    const auto &nni_edge = GetDAG().GetDAGEdge(nni_edge_id);
    for (const auto node_id : {nni_edge.GetParent(), nni_edge.GetChild()}) {
      const auto &node = GetDAG().GetDAGNode(node_id);
      for (const auto dir : DirectionEnum::Iterator()) {
        for (const auto clade : SubsplitCladeEnum::Iterator()) {
          for (const auto adj_node_id : node.GetNeighbors(dir, clade)) {
            const auto parent_id = (dir == Direction::Rootward) ? adj_node_id : node_id;
            const auto child_id = (dir == Direction::Rootward) ? node_id : adj_node_id;
            const auto edge_id = GetDAG().GetEdgeIdx(parent_id, child_id);
            const auto &edge = GetDAG().GetDAGEdge(edge_id);
            SafeAddOutputNNIsToAdjacentNNIs(
                GetDAG().GetDAGNodeBitset(parent_id),
                GetDAG().GetDAGNodeBitset(child_id),
                (edge.GetSubsplitClade() == SubsplitClade::Left));
          }
        }
      }
    }
  }
}

void NNIEngine::UpdateAcceptedNNIs(const bool save_past_nnis) {
  if (save_past_nnis) {
    accepted_past_nnis_.insert(accepted_nnis_.begin(), accepted_nnis_.end());
  }
  if (HasTPEvalEngine()) {
    for (const auto &nni : accepted_nnis_) {
      GetTPEvalEngine().GetScoredNNIs().erase(nni);
    }
  }
  if (HasGPEvalEngine()) {
    for (const auto &nni : accepted_nnis_) {
      GetGPEvalEngine().GetScoredNNIs().erase(nni);
    }
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
  if (GetRescoreRejectedNNIs()) {
    if (HasTPEvalEngine()) {
      GetTPEvalEngine().GetScoredNNIs().clear();
    }
    if (HasGPEvalEngine()) {
      GetGPEvalEngine().GetScoredNNIs().clear();
    }
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
