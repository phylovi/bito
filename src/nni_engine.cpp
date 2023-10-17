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
  const auto it_1 = GetScoredNNIs().find(nni);
  if (it_1 != GetScoredNNIs().end()) {
    return it_1->second;
  }
  const auto it_2 = GetPastScoredNNIs().find(nni);
  if (it_2 != GetScoredNNIs().end()) {
    return it_2->second;
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
  ResetNNIData();
  os << "RunInit::ResetNNIData: " << timer.Lap() << std::endl;
  SyncAdjacentNNIsWithDAG(true);
  os << "RunInit::SyncAdjacentNNIsWithDAG: " << timer.Lap() << std::endl;
  FilterInit();
  os << "RunInit::FilterInit: " << timer.Lap() << std::endl;
}

void NNIEngine::RunMainLoop(const bool is_quiet) {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cout);
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);

  // (1) Add all adjacent NNIs to the GraftDAG.
  GraftAdjacentNNIsToDAG(is_quiet);
  os << "RunMainLoop::GraftAdjacentNNIsToDAG: " << timer.Lap() << std::endl;
  // (2) Compute each adjacent NNI score.
  FilterPreScore();
  os << "RunMainLoop::FilterPreScore: " << timer.Lap() << std::endl;
  // (2b) Optional per-NNI function.
  FilterScoreAdjacentNNIs();
  os << "RunMainLoop::FilterScoreAdjacentNNIs: " << timer.Lap() << std::endl;
  FilterPostScore();
  os << "RunMainLoop::FilterScoreUpdate: " << timer.Lap() << std::endl;
  // (3) Select whether to accept or reject adjacent NNIs via filter.
  FilterEvaluateAdjacentNNIs();
  os << "RunMainLoop::FilterEvaluate: " << timer.Lap() << std::endl;
  // (4a) Remove adjacent NNIs from GraftDAG.
  RemoveAllGraftedNNIsFromDAG();
  // (4b) Add accepted NNIs permanently to DAG.
  AddAcceptedNNIsToDAG(is_quiet);
  os << "RunMainLoop::RemoveAndAddNNIs: " << timer.Lap() << std::endl;

  iter_count_++;
}

void NNIEngine::RunPostLoop(const bool is_quiet) {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cout);
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);

  // (5) Update NNI data.
  // (5a) Update rejected NNIs.
  UpdateRejectedNNIs();
  os << "RunPostLoop::UpdateRejectedNNIs: " << timer.Lap() << std::endl;
  // (5b) Update adjacent and new NNIs.
  UpdateAdjacentNNIs();
  os << "RunPostLoop::UpdateAdjacentNNIs: " << timer.Lap() << std::endl;
  // (5c) Update scored NNIs.
  UpdateScoredNNIs();
  os << "RunPostLoop::UpdateScoredNNIs: " << timer.Lap() << std::endl;
  // (5d) Update accepted NNIs.
  UpdateAcceptedNNIs();
  os << "RunPostLoop::UpdateAcceptedNNIs: " << timer.Lap() << std::endl;
}

// ** Filter Subroutines

void NNIEngine::FilterInit() {
  if (filter_init_fn_) {
    (filter_init_fn_)(*this);
  }
}

void NNIEngine::FilterPreScore() {
  if (filter_pre_score_fn_) {
    (filter_pre_score_fn_)(*this);
  }
}

void NNIEngine::FilterScoreAdjacentNNIs() {
  if (filter_score_loop_fn_) {
    for (const auto &nni : GetNNIsToRescore()) {
      const double nni_score = (filter_score_loop_fn_)(*this, nni);
      AddNNIScore(nni, nni_score);
    }
  }
}

void NNIEngine::FilterPostScore() {
  if (filter_post_score_fn_) {
    (filter_post_score_fn_)(*this);
  }
  Assert(new_sorted_scored_nnis_.size() == GetNNIsToRescore().size(),
         "After scoring NNIs, not all NNIs have been assigned scores.");
}

void NNIEngine::FilterEvaluateAdjacentNNIs() {
  Assert(filter_evaluate_fn_ or filter_evaluate_loop_fn_,
         "Must assign a filter process function before running NNIEngine.");
  if (filter_evaluate_fn_) {
    (filter_evaluate_fn_)(*this, GetNNIsToReevaluate(), GetScoredNNIsToReevaluate(),
                          GetSortedScoredNNIsToReevaluate(), accepted_nnis_);
  }
  if (filter_evaluate_loop_fn_) {
    for (const auto &[nni, nni_score] : GetScoredNNIsToReevaluate()) {
      const bool accept_nni = (filter_evaluate_loop_fn_)(*this, nni, nni_score);
      if (accept_nni) {
        accepted_nnis_.insert(nni);
      } else if (track_rejected_nnis_) {
        rejected_nnis_.insert(nni);
      }
    }
  }
}

void NNIEngine::FilterPostModification(
    const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni) {
  if (filter_post_modification_fn_) {
    (filter_post_modification_fn_)(*this, mods_, nni_to_pre_nni);
  }
}

// ** Filter Setters

void NNIEngine::SetFilterInitFunction(StaticFilterInitFunction filter_init_fn) {
  filter_init_fn_ = filter_init_fn;
}

void NNIEngine::SetFilterPreScoreFunction(
    StaticFilterUpdateFunction filter_pre_update_fn) {
  filter_pre_score_fn_ = filter_pre_update_fn;
}

void NNIEngine::SetFilterScoreLoopFunction(
    StaticFilterScoreLoopFunction filter_score_loop_fn) {
  filter_score_loop_fn_ = filter_score_loop_fn;
}

void NNIEngine::SetFilterPostScoreFunction(
    StaticFilterUpdateFunction filter_post_score_fn) {
  filter_post_score_fn_ = filter_post_score_fn;
}

void NNIEngine::SetFilterEvaluateFunction(
    StaticFilterEvaluateFunction filter_evaluate_fn) {
  filter_evaluate_fn_ = filter_evaluate_fn;
}

void NNIEngine::SetFilterEvaluateLoopFunction(
    StaticFilterEvaluateLoopFunction filter_evaluate_loop_fn) {
  filter_evaluate_loop_fn_ = filter_evaluate_loop_fn;
}

void NNIEngine::SetFilterPostModificationFunction(
    StaticFilterModificationFunction filter_post_modification_fn) {
  filter_post_modification_fn_ = filter_post_modification_fn;
}

// ** Filtering Schemes

void NNIEngine::SetNoFilter(const bool accept_all_nnis) {
  SetNoEvaluate(accept_all_nnis);
}

void NNIEngine::SetGPLikelihoodCutoffFilteringScheme(const double score_cutoff) {
  Assert(HasGPEvalEngine(), "Must MakeGPEvalEngine before using the filtering scheme.");
  SelectGPEvalEngine();
  SetScoreViaEvalEngine();
  SetEvaluateViaMinScoreCutoff(score_cutoff);
}

void NNIEngine::SetTPLikelihoodCutoffFilteringScheme(const double score_cutoff) {
  Assert(HasTPEvalEngine() && GetTPEngine().HasLikelihoodEvalEngine(),
         "Must MakeTPEvalEngine before using the filtering scheme.");
  SelectTPLikelihoodEvalEngine();
  SetScoreViaEvalEngine();
  SetEvaluateViaMinScoreCutoff(score_cutoff);
}

void NNIEngine::SetTPParsimonyCutoffFilteringScheme(const double score_cutoff) {
  Assert(HasTPEvalEngine() && GetTPEngine().HasParsimonyEvalEngine(),
         "Must MakeTPEvalEngine before using the filtering scheme.");
  SelectTPParsimonyEvalEngine();
  SetScoreViaEvalEngine();
  SetEvaluateViaMaxScoreCutoff(score_cutoff);
}

void NNIEngine::SetGPLikelihoodDropFilteringScheme(const double score_cutoff) {
  Assert(HasGPEvalEngine(), "Must MakeGPEvalEngine before using the filtering scheme.");
  SelectGPEvalEngine();
  SetScoreViaEvalEngine();
  SetFilterPostScoreFunction([score_cutoff](NNIEngine &this_nni_engine) {
    double max = this_nni_engine.GetMaxScore();
    this_nni_engine.SetEvaluateViaMinScoreCutoff(max - score_cutoff);
  });
}

void NNIEngine::SetTPLikelihoodDropFilteringScheme(const double score_cutoff) {
  Assert(HasTPEvalEngine() && GetTPEngine().HasLikelihoodEvalEngine(),
         "Must MakeTPEvalEngine before using the filtering scheme.");
  SelectTPLikelihoodEvalEngine();
  SetScoreViaEvalEngine();
  SetFilterPostScoreFunction([score_cutoff](NNIEngine &this_nni_engine) {
    double max = this_nni_engine.GetMaxScore();
    this_nni_engine.SetEvaluateViaMinScoreCutoff(max - score_cutoff);
  });
}

void NNIEngine::SetTPParsimonyDropFilteringScheme(const double score_cutoff) {
  Assert(HasTPEvalEngine() && GetTPEngine().HasParsimonyEvalEngine(),
         "Must MakeTPEvalEngine before using the filtering scheme.");
  SelectTPParsimonyEvalEngine();
  SetScoreViaEvalEngine();
  SetFilterPostScoreFunction([score_cutoff](NNIEngine &this_nni_engine) {
    double min = this_nni_engine.GetMinScore();
    this_nni_engine.SetEvaluateViaMaxScoreCutoff(min + score_cutoff);
  });
}

void NNIEngine::SetTopKScoreFilteringScheme(const size_t k, const bool max_is_best) {
  SetFilterPostScoreFunction([k, max_is_best](NNIEngine &this_nni_engine) {
    auto score_cutoff =
        max_is_best ? this_nni_engine.GetMaxKScore(k) : this_nni_engine.GetMinKScore(k);
    if (max_is_best) {
      this_nni_engine.SetEvaluateViaMinScoreCutoff(score_cutoff);
    } else {
      this_nni_engine.SetEvaluateViaMaxScoreCutoff(score_cutoff);
    }
  });
}

// ** Filtering Scheme Helper Functions

void NNIEngine::SetScoreToConstant(const double value) {
  SetFilterScoreLoopFunction(
      [value](NNIEngine &this_nni_engine, const NNIOperation &nni) { return value; });
}

void NNIEngine::SetScoreViaEvalEngine() {
  SetFilterInitFunction(
      [](NNIEngine &this_nni_engine) { this_nni_engine.PrepEvalEngine(); });
  SetFilterPreScoreFunction(
      [](NNIEngine &this_nni_engine) { this_nni_engine.ScoreAdjacentNNIs(); });
  SetFilterScoreLoopFunction([](NNIEngine &this_nni_engine, const NNIOperation &nni) {
    return this_nni_engine.GetEvalEngine().GetScoredNNIs().find(nni)->second;
  });
  SetFilterPostModificationFunction(
      [](NNIEngine &this_nni_engine, const SubsplitDAG::ModificationResult &mods,
         const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni) {
        this_nni_engine.GrowEvalEngineForDAG(mods.node_reindexer, mods.edge_reindexer);
        this_nni_engine.UpdateEvalEngineAfterModifyingDAG(
            nni_to_pre_nni, mods.prv_node_count, mods.node_reindexer,
            mods.prv_edge_count, mods.edge_reindexer);
      });
}

void NNIEngine::SetNoEvaluate(const bool set_all_nni_to_pass) {
  SetFilterEvaluateFunction(
      [set_all_nni_to_pass](NNIEngine &this_nni_engine, const NNISet &nnis_to_evaluate,
                            const NNIDoubleMap &scored_nnis,
                            const DoubleNNIPairSet &sorted_scored_nnis,
                            NNISet &accepted_nnis) {
        if (set_all_nni_to_pass) {
          for (auto it = sorted_scored_nnis.rbegin(); it != sorted_scored_nnis.rend();
               it++) {
            const auto &[nni_score, nni] = *it;
            accepted_nnis.insert(nni);
          }
        }
      });
}

void NNIEngine::SetEvaluateViaSetOfNNIs(const std::set<NNIOperation> &nnis_to_accept) {
  SetFilterEvaluateFunction([nnis_to_accept](NNIEngine &this_nni_engine,
                                             const NNISet &nnis_to_evaluate,
                                             const NNIDoubleMap &scored_nnis,
                                             const DoubleNNIPairSet &sorted_scored_nnis,
                                             NNISet &accepted_nnis) {
    for (const auto &nni : nnis_to_accept) {
      accepted_nnis.insert(nni);
    }
  });
}

void NNIEngine::SetEvaluateViaMinScoreCutoff(const double score_cutoff) {
  SetFilterEvaluateFunction([score_cutoff](NNIEngine &this_nni_engine,
                                           const NNISet &nnis_to_evaluate,
                                           const NNIDoubleMap &scored_nnis,
                                           const DoubleNNIPairSet &sorted_scored_nnis,
                                           NNISet &accepted_nnis) {
    for (auto it = sorted_scored_nnis.rbegin(); it != sorted_scored_nnis.rend(); it++) {
      const auto &[nni_score, nni] = *it;
      if (nni_score < score_cutoff) {
        break;
      }
      accepted_nnis.insert(nni);
    }
  });
}

void NNIEngine::SetEvaluateViaMaxScoreCutoff(const double score_cutoff) {
  SetFilterEvaluateFunction([score_cutoff](NNIEngine &this_nni_engine,
                                           const NNISet &nnis_to_evaluate,
                                           const NNIDoubleMap &scored_nnis,
                                           const DoubleNNIPairSet &sorted_scored_nnis,
                                           NNISet &accepted_nnis) {
    for (auto it = sorted_scored_nnis.begin(); it != sorted_scored_nnis.end(); it++) {
      const auto &[nni_score, nni] = *it;
      if (nni_score > score_cutoff) {
        break;
      }
      accepted_nnis.insert(nni);
    }
  });
}

void NNIEngine::ScoreAdjacentNNIs() {
  if (IsEvalEngineInUse(NNIEvalEngineType::GPEvalEngine)) {
    GetGPEvalEngine().ScoreAdjacentNNIs(GetNNIsToRescore());
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaLikelihood)) {
    GetTPEvalEngine().ScoreAdjacentNNIs(GetNNIsToRescore());
  }
  if (IsEvalEngineInUse(NNIEvalEngineType::TPEvalEngineViaParsimony)) {
    GetTPEvalEngine().ScoreAdjacentNNIs(GetNNIsToRescore());
  }
}

double NNIEngine::GetMinScore() const {
  const auto &scores = GetSortedScoredNNIsToReevaluate();
  return scores.empty() ? INFINITY : scores.begin()->first;
}

double NNIEngine::GetMaxScore() const {
  const auto &scores = GetSortedScoredNNIsToReevaluate();
  return scores.empty() ? -INFINITY : scores.rbegin()->first;
}

double NNIEngine::GetMinKScore(const size_t k) const {
  const auto &scores = GetSortedScoredNNIsToReevaluate();
  double score = INFINITY;
  if (scores.empty() or k == 0) return score;
  size_t count = 0;
  for (auto it = scores.begin(); it != scores.end(); it++) {
    const auto &[nni_score, nni] = *it;
    score = nni_score;
    if (++count >= k) break;
  }
  return score;
}

double NNIEngine::GetMaxKScore(const size_t k) const {
  const auto &scores = GetSortedScoredNNIsToReevaluate();
  double score = -INFINITY;
  if (scores.empty() or k == 0) return score;
  size_t count = 0;
  for (auto it = scores.rbegin(); it != scores.rend(); it++) {
    const auto &[nni_score, nni] = *it;
    score = nni_score;
    if (++count >= k) break;
  }
  return score;
}

std::set<NNIOperation> NNIEngine::GetMinKScoringNNIs(const size_t k) const {
  const auto &scores = GetSortedScoredNNIsToReevaluate();
  std::set<NNIOperation> nnis;
  if (scores.empty() or k == 0) return nnis;
  size_t count = 0;
  for (auto it = scores.begin(); it != scores.end(); it++) {
    auto &[nni_score, nni] = *it;
    nnis.insert(nni);
    if (++count >= k) break;
  }
  return nnis;
}

std::set<NNIOperation> NNIEngine::GetMaxKScoringNNIs(const size_t k) const {
  const auto &scores = GetSortedScoredNNIsToReevaluate();
  std::set<NNIOperation> nnis;
  if (scores.empty() or k == 0) return nnis;
  size_t count = 0;
  for (auto it = scores.rbegin(); it != scores.rend(); it++) {
    auto &[nni_score, nni] = *it;
    nnis.insert(nni);
    if (++count >= k) break;
  }
  return nnis;
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
  mods_.Reinit(GetDAG());
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
  os << "AddAcceptedNNIsToDAG::Prep: " << timer.Lap() << std::endl;
  // BitsetPairVector nodes_to_add;
  for (const auto &nni : GetAcceptedNNIs()) {
    auto mods = GetDAG().AddNodePair(nni);
    mods_ = mods_.ComposeWith(mods);
    // nodes_to_add.push_back({nni.GetParent(), nni.GetChild()});
  }
  // mods_ = GetDAG().AddNodes(nodes_to_add);
  os << "AddAcceptedNNIsToDAG::AddAllAcceptedNNIs: " << timer.Lap() << std::endl;
  RemoveAllGraftedNNIsFromDAG();
  os << "AddAcceptedNNIsToDAG::RemoveGraftedNNIs: " << timer.Lap() << std::endl;
  FilterPostModification(nni_to_pre_nni);
  os << "AddAcceptedNNIsToDAG::FilterPostModification: " << timer.Lap() << std::endl;
}

void NNIEngine::GraftAdjacentNNIsToDAG(const bool is_quiet) {
  BitsetPairVector nodes_to_add;
  for (const auto &nni : GetNNIsToRescore()) {
    GetGraftDAG().AddNodePair(nni);
    // nodes_to_add.push_back({nni.GetParent(), nni.GetChild()});
  }
  // GetGraftDAG().AddNodes(nodes_to_add);
}

void NNIEngine::RemoveAllGraftedNNIsFromDAG() { GetGraftDAG().RemoveAllGrafts(); }

// ** NNI Maintenance

void NNIEngine::SyncAdjacentNNIsWithDAG(const bool on_init) {
  adjacent_nnis_.clear();
  new_adjacent_nnis_.clear();
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
  // If not on initial run, remove new NNIs that have already been seen.
  if (!on_init) {
    for (const auto &nni : GetPastAcceptedNNIs()) {
      new_adjacent_nnis_.erase(nni);
    }
  }
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
      new_adjacent_nnis_.insert(new_nni);
    }
  }
}

void NNIEngine::AddNNIScore(const NNIOperation &nni, const double score) {
  RemoveNNIScore(nni);
  scored_nnis_.insert({nni, score});
  new_scored_nnis_.insert({nni, score});
  sorted_scored_nnis_.insert({score, nni});
  new_sorted_scored_nnis_.insert({score, nni});
};

void NNIEngine::RemoveNNIScore(const NNIOperation &nni) {
  double old_score = -INFINITY;
  auto it = scored_nnis_.find(nni);
  if (it != scored_nnis_.end()) {
    old_score = it->second;
    sorted_scored_nnis_.erase({old_score, nni});
    new_sorted_scored_nnis_.erase({old_score, nni});
    scored_nnis_.erase(nni);
    new_scored_nnis_.erase(nni);
  }
}

void NNIEngine::UpdateAdjacentNNIs() {
  new_adjacent_nnis_.clear();
  for (const auto &nni : accepted_nnis_) {
    adjacent_nnis_.erase(nni);
    RemoveNNIScore(nni);
  }
  for (const auto &edge_id : mods_.added_edge_idxs) {
    auto nni = GetDAG().GetNNI(edge_id);
    adjacent_nnis_.erase(nni);
    RemoveNNIScore(nni);
  }
  for (const auto &nni : GetAcceptedNNIs()) {
    const auto nni_edge_id = GetDAG().GetEdgeIdx(nni);
    const auto &nni_edge = GetDAG().GetDAGEdge(nni_edge_id);
    for (const auto node_id : {nni_edge.GetParent(), nni_edge.GetChild()}) {
      const auto &node = GetDAG().GetDAGNode(node_id);
      for (const auto dir : DirectionEnum::Iterator()) {
        for (const auto clade : SubsplitCladeEnum::Iterator()) {
          const auto node_view = node.GetNeighbors(dir, clade);
          for (auto it = node_view.begin(); it != node_view.end(); ++it) {
            const auto &edge = GetDAG().GetDAGEdge(it.GetEdgeId());
            SafeAddOutputNNIsToAdjacentNNIs(
                GetDAG().GetDAGNodeBitset(edge.GetParent()),
                GetDAG().GetDAGNodeBitset(edge.GetChild()),
                (edge.GetSubsplitClade() == SubsplitClade::Left));
          }
        }
      }
    }
  }

  if (rescore_old_nnis_adjacent_to_new_nnis_) {
    std::cout << "new_adjacent_nnis [before]: " << new_adjacent_nnis_.size()
              << std::endl;
    UpdateOutOfDateAdjacentNNIs();
    std::cout << "new_adjacent_nnis [after]: " << new_adjacent_nnis_.size()
              << std::endl;
  }
}

void NNIEngine::UpdateOutOfDateAdjacentNNIs() {
  // If there have been modifications surrounding current adjacent NNIs, re-label them
  // as new.
  auto &new_node_ids = mods_.added_node_ids;
  auto &new_edge_ids = mods_.added_edge_idxs;
  std::set<NodeId> updated_node_ids(new_node_ids.begin(), new_node_ids.end());
  for (const auto edge_id : new_edge_ids) {
    const auto edge = GetDAG().GetDAGEdge(edge_id);
    updated_node_ids.insert(edge.GetParent());
    updated_node_ids.insert(edge.GetChild());
  }

  for (const auto &nni : GetNNIsToReevaluate()) {
    // If NNI is already new, skip.
    if (new_adjacent_nnis_.find(nni) != new_adjacent_nnis_.end()) {
      continue;
    }
    bool nni_found = false;
    for (const auto &nni_subsplit : {nni.GetParent(), nni.GetChild()}) {
      if (nni_found) {
        continue;
      }
      // First, if one of the nodes in the NNI are newly added, re-label as new.
      if (GetDAG().ContainsNode(nni_subsplit)) {
        const auto node_id = GetDAG().GetDAGNodeId(nni_subsplit);
        if (updated_node_ids.find(node_id) != updated_node_ids.end()) {
          new_adjacent_nnis_.insert(nni);
          continue;
        }
      }
      // Second, if any edges connecting to grandparent or grandchild nodes are new,
      // re-label as new.
      const auto [left_parent_ids, right_parent_ids] =
          GetDAG().FindParentNodeIdsViaMap(nni_subsplit);
      const auto [left_child_ids, right_child_ids] =
          GetDAG().FindChildNodeIdsViaMap(nni_subsplit);
      for (const auto &node_ids :
           {left_parent_ids, right_parent_ids, left_child_ids, right_child_ids}) {
        if (nni_found) break;
        for (const auto node_id : node_ids) {
          if (updated_node_ids.find(node_id) != updated_node_ids.end()) {
            new_adjacent_nnis_.insert(nni);
            nni_found = true;
            break;
          }
        }
      }
    }
  }
}

void NNIEngine::UpdateRejectedNNIs() {
  if (save_past_rejected_nnis_) {
    rejected_past_nnis_.insert(new_adjacent_nnis_.begin(), new_adjacent_nnis_.end());
    for (const auto &nni : accepted_nnis_) {
      rejected_past_nnis_.erase(nni);
    }
  }
  rejected_nnis_.clear();
}

void NNIEngine::UpdateScoredNNIs() {
  if (save_past_scored_nnis_) {
    scored_past_nnis_.insert(scored_nnis_.begin(), scored_nnis_.end());
    for (const auto &nni : accepted_nnis_) {
      scored_past_nnis_.erase(nni);
    }
  }
  new_scored_nnis_.clear();
  new_sorted_scored_nnis_.clear();
}

void NNIEngine::UpdateAcceptedNNIs() {
  if (save_past_accepted_nnis_) {
    accepted_past_nnis_.insert(accepted_nnis_.begin(), accepted_nnis_.end());
  }
  accepted_nnis_.clear();
}

void NNIEngine::ResetNNIData() {
  new_adjacent_nnis_.clear();
  adjacent_nnis_.clear();
  accepted_nnis_.clear();
  accepted_past_nnis_.clear();
  rejected_nnis_.clear();
  rejected_past_nnis_.clear();
}
