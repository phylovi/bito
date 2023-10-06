// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// NNI Engine is used to explore the topological space surrounding the
// subsplitDAG using NNI (Nearest Neighbor Interchange) moves. The engine functions by
// finding all the NNIs adjacent to the DAG, scoring them, then choosing to accept or
// reject each NNI based on a filtering criteria.  For each NNI added to the DAG, each
// new edge added also adds new adjacent NNIs to the DAG.  This process is repeated
// until no NNIs pass the filtering criteria. NNI scoring is facilitated by the
// NNIEvalEngines.
//
// To understand this code, it may be easiest to start by reading RunMainLoop and then
// branching out from there.

#pragma once

#include "gp_engine.hpp"
#include "tp_engine.hpp"
#include "gp_dag.hpp"

#include "bitset.hpp"
#include "subsplit_dag.hpp"
#include "nni_operation.hpp"
#include "graft_dag.hpp"
#include "sugar.hpp"
#include "gp_operation.hpp"
#include "reindexer.hpp"
#include "nni_evaluation_engine.hpp"
#include "nni_engine_key_index.hpp"

enum class NNIEvalEngineType {
  GPEvalEngine,
  TPEvalEngineViaLikelihood,
  TPEvalEngineViaParsimony
};
static const inline size_t NNIEvalEngineTypeCount = 3;
class NNIEvalEngineTypeEnum
    : public EnumWrapper<NNIEvalEngineType, size_t, NNIEvalEngineTypeCount,
                         NNIEvalEngineType::GPEvalEngine,
                         NNIEvalEngineType::TPEvalEngineViaParsimony> {
 public:
  static inline const std::string Prefix = "EvalEngineType";
  static inline const Array<std::string> Labels = {{"GPEvalEngine", "TPEvalEngine"}};

  static std::string ToString(const NNIEvalEngineType e) {
    std::stringstream ss;
    ss << Prefix << "::" << Labels[e];
    return ss.str();
  }
  friend std::ostream &operator<<(std::ostream &os, const NNIEvalEngineType e) {
    os << ToString(e);
    return os;
  }
};

class NNIEngine {
 public:
  using NNIDoubleMap = std::map<NNIOperation, double>;
  using DoubleNNIPairSet = std::set<std::pair<double, NNIOperation>>;

  // Constructors
  NNIEngine(GPDAG &dag, std::optional<GPEngine *> gp_engine = std::nullopt,
            std::optional<TPEngine *> tp_engine = std::nullopt);

  // ** Access

  // Get Reference of DAG.
  GPDAG &GetDAG() { return dag_; };
  const GPDAG &GetDAG() const { return dag_; };
  // Get Reference of GraftDAG.
  GraftDAG &GetGraftDAG() { return *graft_dag_.get(); };
  const GraftDAG &GetGraftDAG() const { return *graft_dag_.get(); };
  // Get Reference of Evaluation Engine.
  NNIEvalEngine &GetEvalEngine() {
    Assert(HasEvalEngine(), "EvalEngine has not been set.");
    return *eval_engine_;
  }
  const NNIEvalEngine &GetEvalEngine() const {
    Assert(HasEvalEngine(), "EvalEngine has not been set.");
    return *eval_engine_;
  }
  bool HasEvalEngine() const { return eval_engine_ != nullptr; }
  // Get Reference of GPEvalEngine.
  NNIEvalEngineViaGP &GetGPEvalEngine() {
    Assert(HasGPEvalEngine(), "GPEvalEngine has not been set.");
    return *eval_engine_via_gp_.get();
  }
  const NNIEvalEngineViaGP &GetGPEvalEngine() const {
    Assert(HasGPEvalEngine(), "GPEvalEngine has not been set.");
    return *eval_engine_via_gp_.get();
  }
  bool HasGPEvalEngine() const { return eval_engine_via_gp_ != nullptr; }
  // Get Reference of TPEvalEngine.
  NNIEvalEngineViaTP &GetTPEvalEngine() {
    Assert(HasTPEvalEngine(), "TPEvalEngine has not been set.");
    return *eval_engine_via_tp_.get();
  }
  const NNIEvalEngineViaTP &GetTPEvalEngine() const {
    Assert(HasTPEvalEngine(), "TPEvalEngine has not been set.");
    return *eval_engine_via_tp_.get();
  }
  bool HasTPEvalEngine() const { return eval_engine_via_tp_ != nullptr; }
  // Get Reference of GPEngine.
  const GPEngine &GetGPEngine() const {
    Assert(HasGPEvalEngine(), "GPEvalEngine has not been set.");
    return GetGPEvalEngine().GetGPEngine();
  }
  // Get Reference of TPEngine.
  const TPEngine &GetTPEngine() const {
    Assert(HasTPEvalEngine(), "TPEvalEngine has not been set.");
    return GetTPEvalEngine().GetTPEngine();
  }

  // Get score for given NNI in or adjacent to DAG.
  double GetScoreByNNI(const NNIOperation &nni) const;
  double GetScoreByEdge(const EdgeId edge_id) const;

  // NNIs currently adjacent to DAG.
  const NNISet &GetAdjacentNNIs() const { return adjacent_nnis_; }
  size_t GetAdjacentNNICount() const { return adjacent_nnis_.size(); }
  // Adjacent NNIs that have been added on the current iteration.
  const NNISet &GetNewAdjacentNNIs() const { return new_adjacent_nnis_; }
  size_t GetNewAdjacentNNICount() const { return new_adjacent_nnis_.size(); }
  size_t GetOldNNICount() const {
    return adjacent_nnis_.size() - new_adjacent_nnis_.size();
  }
  // NNIs that have been Accepted on current iteration.
  const NNISet &GetAcceptedNNIs() const { return accepted_nnis_; }
  size_t GetAcceptedNNICount() const { return GetAcceptedNNIs().size(); }
  // NNIs that have been Accepted from all iterations.
  const NNISet &GetPastAcceptedNNIs() const { return accepted_past_nnis_; }
  size_t GetPastAcceptedNNICount() const { return GetPastAcceptedNNIs().size(); }
  // Get NNIs that have been Rejected on current iteration.
  const NNISet &GetRejectedNNIs() const { return rejected_nnis_; }
  size_t GetRejectedNNICount() const { return GetRejectedNNIs().size(); }
  // Get NNIs that have been Rejected from all iterations.
  const NNISet &GetPastRejectedNNIs() const { return rejected_past_nnis_; }
  size_t GetPastRejectedNNICount() const { return GetPastRejectedNNIs().size(); }
  // Get Map of proposed NNIs with their score.
  const NNIDoubleMap &GetScoredNNIs() const { return scored_nnis_; }
  size_t GetScoredNNICount() const { return GetScoredNNIs().size(); }
  // Get Map of proposed NNIs with their score from all iterations.
  const NNIDoubleMap &GetPastScoredNNIs() const { return scored_past_nnis_; }
  size_t GetPastScoredNNICount() const { return GetPastScoredNNIs().size(); }
  // Get NNIs to rescore or NNIs to re-evaluate.
  const NNISet &GetNNIsToRescore() const {
    return GetRescoreRejectedNNIs() ? GetAdjacentNNIs() : GetNewAdjacentNNIs();
  }
  const NNISet &GetNNIsToReevaluate() const {
    return GetReevaluateRejectedNNIs() ? GetAdjacentNNIs() : GetNewAdjacentNNIs();
  }
  const NNIDoubleMap &GetScoredNNIsToRescore() const {
    return GetRescoreRejectedNNIs() ? scored_nnis_ : new_scored_nnis_;
  }
  const NNIDoubleMap &GetScoredNNIsToReevaluate() const {
    return GetReevaluateRejectedNNIs() ? scored_nnis_ : new_scored_nnis_;
  }
  const DoubleNNIPairSet &GetSortedScoredNNIsToRescore() const {
    return GetRescoreRejectedNNIs() ? sorted_scored_nnis_ : new_sorted_scored_nnis_;
  }
  const DoubleNNIPairSet &GetSortedScoredNNIsToReevaluate() const {
    return GetReevaluateRejectedNNIs() ? sorted_scored_nnis_ : new_sorted_scored_nnis_;
  }
  // Get vector of proposed NNI scores.
  DoubleVector GetNNIScores() const {
    DoubleVector scores;
    for (const auto &[nni, score] : GetScoredNNIs()) {
      std::ignore = nni;
      scores.push_back(score);
    }
    return scores;
  }
  // Get proposed NNI score.
  double GetNNIScore(const NNIOperation &nni) const {
    const auto it_1 = GetScoredNNIs().find(nni);
    if (it_1 != GetScoredNNIs().end()) {
      return it_1->second;
    }
    const auto it_2 = GetPastScoredNNIs().find(nni);
    if (it_2 != GetPastScoredNNIs().end()) {
      return it_2->second;
    }
    // Failwith("ERROR: Could not find NNI in ScoredNNIs.");
    return -INFINITY;
  }

  // Reindexers for recent DAG modifications.
  const Reindexer &GetNodeReindexer() const { return mods_.node_reindexer; }
  const Reindexer &GetEdgeReindexer() const { return mods_.edge_reindexer; }

  // Option whether to re-evaluate rejected nnis.
  bool GetReevaluateRejectedNNIs() const { return reevaluate_rejected_nnis_; }
  void SetReevaluateRejectedNNIs(const bool reevaluate_rejected_nnis) {
    reevaluate_rejected_nnis_ = reevaluate_rejected_nnis;
  }
  // Option whether to re-score rejected nnis.
  bool GetRescoreRejectedNNIs() const { return rescore_rejected_nnis_; }
  void SetRescoreRejectedNNIs(const bool rescore_rejected_nnis) {
    rescore_rejected_nnis_ = rescore_rejected_nnis;
  }
  // Option whether to include NNIs at containing rootsplits.
  bool GetIncludeRootsplitNNIs() const { return include_rootsplit_nnis_; }
  void SetIncludeRootsplitNNIs(const bool include_rootsplit_nnis) {
    include_rootsplit_nnis_ = include_rootsplit_nnis;
  }

  // Get number of runs of NNI engine.
  size_t GetIterationCount() const { return iter_count_; };
  // Reset number of iterations.
  void ResetIterationCount() { iter_count_ = 0; }

  // ** NNI Evaluation Engine

  // Set GP Engine.
  NNIEvalEngineViaGP &MakeGPEvalEngine(GPEngine *gp_engine);
  // Set TP Engine.
  NNIEvalEngineViaTP &MakeTPEvalEngine(TPEngine *tp_engine);
  // Check if evaluation engine is currently in use.
  bool IsEvalEngineInUse(const NNIEvalEngineType eval_engine_type) {
    return eval_engine_in_use_[eval_engine_type];
  }
  // Remove all evaluation engines from use.
  void ClearEvalEngineInUse();
  // Set evaluation engine type for use in runner.
  void SelectEvalEngine(const NNIEvalEngineType eval_engine_type);
  // Set GP evaluation engine for use in runner.
  void SelectGPEvalEngine();
  // Set TP likelihood evaluation engine for use in runner.
  void SelectTPLikelihoodEvalEngine();
  // Set TP parsimony evaluation engine for use in runner.
  void SelectTPParsimonyEvalEngine();

  // Initial GPEngine for use with GraftDAG.
  void InitEvalEngine();
  // Populate PLVs for quick lookup of likelihoods.
  void PrepEvalEngine();
  // Resize Engine for modified DAG.
  void GrowEvalEngineForDAG(std::optional<Reindexer> node_reindexer,
                            std::optional<Reindexer> edge_reindexer);
  // Update PVs after modifying the DAG.
  void UpdateEvalEngineAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &pre_nni_to_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer);
  // Fetches Pre-NNI data to prep Post-NNI for likelihood computation. Method stores
  // intermediate values in the GPEngine temp space (expects GPEngine has already been
  // resized).
  void GrowEvalEngineForAdjacentNNIs(const bool via_reference = true,
                                     const bool use_unique_temps = false);

  // ** Runners
  // These start the engine, which procedurally ranks and adds (and maybe removes) NNIs
  // to the DAG, until some termination criteria has been satisfied.

  // Primary Runner for NNI Engine.
  void Run(const bool is_quiet = true);
  // Initialization step run before loop.
  void RunInit(const bool is_quiet = true);
  // Step that finds adjacent NNIs, evaluates, then accepts or rejects them.
  void RunMainLoop(const bool is_quiet = true);
  // Step that runs at the end of each loop, preps for next loop.
  void RunPostLoop(const bool is_quiet = true);

  // ** Filter Functions

  // Initialization step for filter before beginning
  using StaticFilterInitFunction = std::function<void(NNIEngine &)>;
  // Update step for any filter computation work to be done at start of each sweep.
  using StaticFilterUpdateFunction = std::function<void(NNIEngine &)>;
  // Function template for computational evaluation to be performed on an adjacent NNI.
  using StaticFilterScoreFunction =
      std::function<double(NNIEngine &, const NNIOperation &)>;
  // Function template for processing an adjacent NNI to be accepted or rejected.
  using StaticFilterEvaluateFunction =
      std::function<void(NNIEngine &, const NNISet &, const NNIDoubleMap &,
                         const DoubleNNIPairSet &, NNISet &)>;
  using StaticFilterEvaluateLoopFunction =
      std::function<bool(NNIEngine &, const NNIOperation &, const double)>;

  // ** Filter Subroutines

  // Initialize filter before first NNI sweep.
  void FilterInit();
  // Update filter parameters for each NNI sweep (before evaluation).
  void FilterPreUpdate();
  // Perform computation on each Adjacent NNI.
  void FilterScoreAdjacentNNIs();
  // Update filter parameters for each NNI sweep (after evaluation).
  void FilterPostUpdate();
  // Apply the filtering method to determine whether each Adjacent NNI will be added to
  // Accepted NNI or Rejected NNI.
  void FilterEvaluateAdjacentNNIs();

  // Set filter initialization function. Called at the beginning of NNI engine run,
  // before main loop.
  void SetFilterInitFunction(StaticFilterInitFunction filter_init_fn);
  // Set filter pre-eval update step.  Performed once each loop before the filter
  // evaluvation step.
  void SetFilterPreUpdateFunction(StaticFilterUpdateFunction filter_pre_update_fn);
  // Set filter evaluation step.  Evaluation step is performed on each proposed adjacent
  // NNI individually, and returns a score for that NNI.
  void SetFilterScoreFunction(StaticFilterScoreFunction filter_score_fn);
  // Set filter post-eval update step. Performed once each loop after the filter
  // evaluvation step.
  void SetFilterPostUpdateFunction(StaticFilterUpdateFunction filter_post_update_fn);
  // Set filter processing step.  Processing step is performed on each proposed adjacent
  // NNI individually, taking in its NNI score and outputting a boolean whether to
  // accept or reject the NNI.
  void SetFilterEvaluateFunction(StaticFilterEvaluateFunction filter_evaluate_fn);
  void SetFilterEvaluateLoopFunction(
      StaticFilterEvaluateLoopFunction filter_evaluate_loop_fn);

  // ** Filtering Schemes

  // Set filtering scheme to simply accept or reject all NNIs.
  void SetNoFilter(const bool accept_all_nnis);

  // Set filtering scheme to use GP likelihoods, using static cutoff.
  void SetGPLikelihoodCutoffFilteringScheme(const double score_cutoff);
  // Set filtering scheme to use TP likelihoods, using static cutoff.
  void SetTPLikelihoodCutoffFilteringScheme(const double score_cutoff);
  // Set filtering scheme to use TP parsimony, using static cutoff.
  void SetTPParsimonyCutoffFilteringScheme(const double score_cutoff);

  // Set filtering scheme to use GP likelihoods, using static cutoff.
  void SetGPLikelihoodDropFilteringScheme(const double score_cutoff);
  // Set filtering scheme to use TP likelihoods, using static cutoff.
  void SetTPLikelihoodDropFilteringScheme(const double score_cutoff);
  // Set filtering scheme to use TP parsimony, using static cutoff.
  void SetTPParsimonyDropFilteringScheme(const double score_cutoff);

  // Set filtering scheme to find the top K best-scoring NNIs.
  void SetTopKScoreFilteringScheme(const size_t k, const bool max_is_best = true);

  // ** Filtering Scheme Helper Functions

  // Set filter to score to constant value.
  void SetScoreToConstant(const double value = -INFINITY);
  // Set filter to score using evaluation engine.
  void SetScoreViaEvalEngine();
  // Set filter to accept/deny all adjacent NNIs.
  void SetNoEvaluate(const bool set_all_nni_to_pass = true);
  // Set filter by accepting NNIs contained in explicit set.
  void SetEvaluateViaSetOfNNIs(const std::set<NNIOperation> &nnis_to_accept);
  // Set cutoff filter to constant cutoff. Accept scores above threshold.
  void SetEvaluateViaMinScoreCutoff(const double score_cutoff);
  // Set cutoff filter to constant cutoff. Accept scores below threshold.
  void SetEvaluateViaMaxScoreCutoff(const double score_cutoff);

  // Performs entire scoring computation for all Adjacent NNIs.
  void ScoreAdjacentNNIs();
  // Get minimum score from scored NNIs.
  double GetMinScore() const;
  // Get maximum score from scored NNIs.
  double GetMaxScore() const;
  // Get bottom kth score from scored NNIs.
  double GetMinKScore(const size_t k) const;
  // Get top kth score from scored NNIs.
  double GetMaxKScore(const size_t k) const;
  // Get bottom kth score from scored NNIs.
  std::set<NNIOperation> GetMinKScoringNNIs(const size_t k) const;
  // Get top kth score from scored NNIs.
  std::set<NNIOperation> GetMaxKScoringNNIs(const size_t k) const;

  // ** Key Indexing

  using KeyIndex = NNIEngineKeyIndex;
  using KeyIndexPairArray = NNIEngineKeyIndexPairArray;
  using KeyIndexMap = NNIEngineKeyIndexMap;
  using KeyIndexMapPair = NNIEngineKeyIndexMapPair;

  // Translate NNIClade type to corresponding PHat PLV from KeyIndex type.
  using NNIClade = NNIOperation::NNIClade;
  static KeyIndex NNICladeToPHatPLV(NNIClade clade_type);

  // Builds an array containing a mapping of plvs from Pre-NNI to Post-NNI, according to
  // remapping of the NNI clades (sister, right_child, left_child).
  static KeyIndexPairArray BuildKeyIndexTypePairsFromPreNNIToPostNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni);

  // Create map of key indices for given NNIOperation needed for computing NNI
  // Likelihoods.
  template <typename DAGType>
  static KeyIndexMap BuildKeyIndexMapForNNI(const NNIOperation &nni, const DAGType &dag,
                                            const size_t node_count);
  KeyIndexMap BuildKeyIndexMapForNNI(const NNIOperation &nni,
                                     const size_t node_count) const;

  // Create map of key indices for Post-NNI, using Pre-NNI map as a reference.
  template <typename DAGType>
  static KeyIndexMap BuildKeyIndexMapForPostNNIViaReferencePreNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni,
      const KeyIndexMap &pre_key_idx, const DAGType &dag);
  KeyIndexMap BuildKeyIndexMapForPostNNIViaReferencePreNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni,
      const KeyIndexMap &pre_key_idx) const;

  // ** DAG Maintenance

  // Add all Accepted NNIs to Main DAG.
  void AddAcceptedNNIsToDAG(const bool is_quiet = true);
  // Add all Adjacent NNIs to Graft DAG.
  void GraftAdjacentNNIsToDAG(const bool is_quiet = true);
  // Remove all NNIs from Graft DAG.
  void RemoveAllGraftedNNIsFromDAG();

  // ** NNI Maintenance
  // These maintain NNIs to stay consistent with the state of associated GraftDAG.

  // Add score to given NNI.
  void AddNNIScore(const NNIOperation &nni, const double score);
  // Remove score for given NNI.
  void RemoveNNIScore(const NNIOperation &nni);

  // Freshly synchonizes NNISet to match the current state of its DAG. Wipes old NNI
  // data and finds all all parent/child pairs adjacent to DAG by iterating over all
  // internal edges in the DAG. (For each internal edges, two NNIs are possible.)
  void SyncAdjacentNNIsWithDAG(const bool on_init = false);
  // Updates NNI Set after given parent/child node pair have been added to the DAG.
  // Removes pair from NNI Set and adds adjacent pairs coming from newly created edges.
  void UpdateAdjacentNNIsAfterDAGAddNodePair(const NNIOperation &nni);
  void UpdateAdjacentNNIsAfterDAGAddNodePair(const Bitset &parent_bitset,
                                             const Bitset &child_bitset);
  // Adds all NNIs from all (node_id, other_id) pairs, where other_id's are elements of
  // the adjacent_node_ids vector. is_edge_leafward tells whether node_id is the child
  // or parent. is_edge_edgeclade determines which side of parent the child descends
  // from.
  void AddAllNNIsFromNodeVectorToAdjacentNNIs(const NodeId node_id,
                                              const SizeVector &adjacent_node_ids,
                                              const bool is_edge_on_left,
                                              const bool is_edge_leafward);
  // Based on given input NNIOperation, produces the two possible output NNIOperations
  // and adds those results to the NNI Set (if results are not a member of the DAG or
  // NNI Set).
  void SafeAddOutputNNIsToAdjacentNNIs(const Bitset &parent_bitset,
                                       const Bitset &child_bitset,
                                       const bool is_edge_on_left);

  // This handles updating all NNI data: adjacent, new, accepted, rejected, and scored
  // NNIs.
  void UpdateRejectedNNIs();
  void UpdateAdjacentNNIs();
  void UpdateScoredNNIs();
  void UpdateAcceptedNNIs();

  // Reset all NNIs, current and past.
  void ResetNNIData();

 private:
  // ** Access

  // Get Reference of GP Engine.
  GPEngine &GetGPEngine() {
    Assert(HasGPEvalEngine(), "GPEvalEngine has not been set.");
    return GetGPEvalEngine().GetGPEngine();
  }
  // Get Reference of TP Engine.
  TPEngine &GetTPEngine() {
    Assert(HasTPEvalEngine(), "TPEvalEngine has not been set.");
    return GetTPEvalEngine().GetTPEngine();
  }

  // Un-owned reference DAG.
  GPDAG &dag_;
  // For adding temporary NNIs to DAG.
  std::unique_ptr<GraftDAG> graft_dag_;
  // Tracks modifications to the DAG.
  SubsplitDAG::ModificationResult mods_;

  // A map showing which Evaluation Engines are "in use".  Several engines may be
  // instatiated, but may or may not be currently used for computation, and therefore
  // may not need to be upkept.
  NNIEvalEngineTypeEnum::Array<bool> eval_engine_in_use_;
  // Un-owned reference to NNI Evaluation Engine. Can be used to evaluate NNIs according
  // to Generalized Pruning, Likelihood, Parsimony, etc. Primary eval engine reference,
  // points to one of the available evaluation engines.
  NNIEvalEngine *eval_engine_ = nullptr;
  std::vector<NNIEvalEngine *> available_eval_engines_;
  // Evaluation engine for scoring NNIs using GP.
  std::unique_ptr<NNIEvalEngineViaGP> eval_engine_via_gp_ = nullptr;
  // Evaluation engine for scoring NNIs using TP.
  std::unique_ptr<NNIEvalEngineViaTP> eval_engine_via_tp_ = nullptr;

  // Set of NNIs to be evaluated, which are a single NNI.
  NNISet adjacent_nnis_;
  // Set of NNIs new to the current iteration.
  NNISet new_adjacent_nnis_;
  // NNIs which have passed the filtering threshold during current iteration, to be
  // added to the DAG.
  NNISet accepted_nnis_;
  // NNIs which have been accepted in a previous iteration of the search.
  NNISet accepted_past_nnis_;
  // NNIs which have failed the filtering threshold during current iteration, NOT to be
  // added to the DAG.
  NNISet rejected_nnis_;
  // NNIs which have been rejected in a previous iteration of the search.
  NNISet rejected_past_nnis_;

  // Map of adjacent NNIs to their score.
  NNIDoubleMap scored_nnis_;
  // Map of new NNIs to their score.
  NNIDoubleMap new_scored_nnis_;
  // Map of previous rejected NNIs to their score.
  NNIDoubleMap scored_past_nnis_;
  // Set of adjacent NNI scores, sorted by score.
  DoubleNNIPairSet sorted_scored_nnis_;
  // Set of new NNI scores, sorted by score.
  DoubleNNIPairSet new_sorted_scored_nnis_;

  // Steps of filtering scheme.
  StaticFilterInitFunction filter_init_fn_ = nullptr;
  StaticFilterUpdateFunction filter_pre_update_fn_ = nullptr;
  StaticFilterScoreFunction filter_score_fn_ = nullptr;
  StaticFilterUpdateFunction filter_post_update_fn_ = nullptr;
  StaticFilterEvaluateFunction filter_evaluate_fn_ = nullptr;
  StaticFilterEvaluateLoopFunction filter_evaluate_loop_fn_ = nullptr;

  // Count number of loops executed by engine.
  size_t iter_count_ = 0;
  // Count number of proposed NNIs computed.
  size_t proposed_nnis_computed_ = 0;

  // Whether to re-evaluate rejected NNIs from previous iterations.
  bool reevaluate_rejected_nnis_ = true;
  // Whether to re-compute scores for rejected NNIs from previous iterations.
  bool rescore_rejected_nnis_ = false;
  // Whether to include NNIs whose parent is a rootsplit.
  bool include_rootsplit_nnis_ = true;
  // Whether to save past iteration data.
  bool save_past_scored_nnis_ = false;
  bool save_past_accepted_nnis_ = true;
  bool save_past_rejected_nnis_ = true;
  bool track_rejected_nnis_ = false;
};
