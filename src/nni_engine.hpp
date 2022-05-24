// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// NNI Engine
// This engine is used to explore the topological space surrounding the
// subsplitDAG using the smallest topological change of the NNI (Nearest Neighbor
// Interchange). The engine functions by finding all the NNIs adjacent to the DAG,
// evaluating them by their likelihoods, then choosing to accept or reject each NNI
// based on a filtering criteria.  For each NNI added to the DAG, each new edge added
// also adds new adjacent NNIs to the DAG.  This process is repeated until no NNIs pass
// the filtering criteria.

#pragma once

#include "gp_engine.hpp"
#include "gp_dag.hpp"

#include "bitset.hpp"
#include "subsplit_dag.hpp"
#include "nni_operation.hpp"
#include "graft_dag.hpp"
#include "sugar.hpp"
#include "gp_operation.hpp"
#include "reindexer.hpp"

using NNIDoubleMap = std::map<NNIOperation, double>;

class NNIEngine {
 public:
  // Constructors:
  NNIEngine(GPDAG &dag, GPEngine &gp_engine);

  // ** Access

  // Get Reference of DAG.
  const GPDAG &GetGPDAG() const { return dag_; }
  // Get Reference of GraftDAG.
  GraftDAG &GetGraftDAG() const { return *graft_dag_.get(); }
  // Get Reference of GPEngine.
  const GPEngine &GetGPEngine() const { return gp_engine_; }
  // Get Adjacent NNIs to DAG.
  const NNISet &GetAdjacentNNIs() const { return adjacent_nnis_; }
  // Get Number of Adjacent NNIs.
  size_t GetAdjacentNNICount() const { return adjacent_nnis_.size(); }
  // Get NNIs that have been accepted into the DAG on current pass.
  const NNISet &GetAcceptedNNIs() const { return accepted_nnis_; }
  // Get Number of Accepted NNIs.
  size_t GetAcceptedNNICount() const { return accepted_nnis_.size(); }
  // Get all past NNIs that have been accepted into the DAG during run.
  const NNISet &GetPastAcceptedNNIs() const { return accepted_past_nnis_; }
  // Get Number of Past Accepted NNIs.
  size_t GetPastAcceptedNNICount() const { return accepted_past_nnis_.size(); }
  // Get NNIs that have been rejected from the DAG on current pass.
  const NNISet &GetRejectedNNIs() const { return rejected_nnis_; }
  // Get Number of Rejected NNIs.
  size_t GetRejectedNNICount() const { return rejected_nnis_.size(); }
  // Get all past NNIs that have been rejected from the DAG during run.
  const NNISet &GetPastRejectedNNIs() const { return rejected_past_nnis_; }
  // Get Number of Past Rejected NNIs.
  size_t GetPastRejectedNNICount() const { return rejected_past_nnis_.size(); }
  // Get Map of proposed NNIs with their score.
  const NNIDoubleMap &GetScoredNNIs() const { return scored_nnis_; }
  // Get number of iterations of NNI engine.
  size_t GetSweepCount() const { return sweep_count_; }
  // Get DAG likelihoods.
  const EigenVectorXd &GetDAGLikelihoods() const { return *dag_likelihoods_; }
  // Get DAG likelihoods for
  double GetDAGEdgeLikelihood(const size_t edge_idx) const {
    const auto &dag_likelihoods = GetDAGLikelihoods();
    Assert(edge_idx < size_t(dag_likelihoods.size()),
           "Requested likelihood for edge index is out-of-range.");
    return dag_likelihoods[edge_idx];
  }

  // ** Setters

  // Set whether to track accepted past NNIs.
  void SetKeepAcceptedNNIs(const bool keep_accepted_past_nnis) {
    keep_accepted_past_nnis_ = keep_accepted_past_nnis;
  }
  // Set whether to track rejected past NNIs.
  void SetKeepRejectedNNIs(const bool keep_rejected_past_nnis) {
    keep_rejected_past_nnis_ = keep_rejected_past_nnis;
  }
  // Set whether to re-evalute rejected past NNIs.
  void SetReevaluateRejecteddNNIs(const bool reevaluate_rejected_nnis) {
    reevaluate_rejected_nnis_ = reevaluate_rejected_nnis;
  }

  // ** Runners
  // These start the engine, which procedurally ranks and adds (and maybe removes) NNIs
  // to the DAG, until some termination criteria has been satisfied.

  // Primary Runner using GraftDAG.
  void Run();
  // Runner subroutines.
  void RunInit();
  // Graft Adjacent NNIs to GraftDAG. Evaluate and Accept or Reject Adjacent NNIs.
  void RunMainLoop();
  // Remove Graft and add Accepted NNIs to HostDAG. Update Adjacent NNIs to account for
  // new Accepted NNIs.  Update Accepted and Rejected NNI lists.
  void RunPostLoop();

  // ** Runner Subroutines

  // Initialize filter before first NNI sweep.
  void FilterInit();
  // Update filter parameters for each NNI sweep (before evaluation).
  void FilterPreUpdate();
  // Apply the filtering method to determine whether each Adjacent NNI will be added
  // to Accepted NNI or Rejected NNI.
  void FilterProcessAdjacentNNIs();
  // Update filter parameters for each NNI sweep (after evaluation).
  void FilterPostUpdate();

  // ** Static Filter Functions

  // Initialization step for filter before beginning engine run.
  using StaticFilterInitFunction =
      std::function<void(NNIEngine &, GPEngine &, GraftDAG &)>;
  // Update step for any filter computation work to be done at start of each sweep.
  using StaticFilterUpdateFunction =
      std::function<void(NNIEngine &, GPEngine &, GraftDAG &)>;
  // Function template for processing an adjacent NNI to be accepted or rejected.
  using StaticFilterProcessFunction =
      std::function<bool(NNIEngine &, GPEngine &, GraftDAG &, const NNIOperation &)>;

  // This StaticFilterInitFunction does nothing.
  static void NoInit(NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                     GraftDAG &this_graft_dag);
  // This StaticFilterUpdateFunction does nothing.
  static void NoUpdate(NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                       GraftDAG &this_graft_dag);
  // This StaticFilterProcessFunction accepts all or rejects all NNIs, according to
  // templated parameter.
  template <bool accept = true>
  static bool NoFilter(NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                       GraftDAG &this_graft_dag, const NNIOperation &nni) {
    return accept;
  }
  // This StaticFilterProcessFunction is for Hill-climbing, only accepting a proposed
  // NNI if it has an improved likelihood over Pre-NNI.
  static bool FilterPostNNIBetterThanPreNNI(NNIEngine &this_nni_engine,
                                            GPEngine &this_gp_engine,
                                            GraftDAG &this_graft_dag,
                                            const NNIOperation &nni);
  // This StaticFilterUpdateFunction computes likelihood for all proposed NNIs.
  static void UpdateScoreAdjacentNNIsByLikelihood(NNIEngine &this_nni_engine,
                                                  GPEngine &this_gp_engine,
                                                  GraftDAG &this_graft_dag);

  // Finds best Pre-NNI likelihood in DAG for given NNI in HostDAG.  For use in
  // hill-climbing filtering scheme.
  // - NOTE: For any proposed NNI in the GraftDAG, there are two possible NNIs which
  // could have produced it, depending on which child clade is swapped with the parent
  // sister clade.  One or both of those Pre-NNIs could be present in the HostDAG.  This
  // function checks for both, and if both are present, selects the one with the highest
  // likelihood.
  double GetBestPreNNILikelihood(const NNIOperation &nni);
  // Get Post-NNI likelihood for proposed NNI in GraftDAG.
  double GetPostNNILikelihood(const NNIOperation &nni);

  // ** NNI Likelihoods
  // Functions for computation or estimating likelihood of NNIs.

  // Map for storing important key indices for NNI Likelihood computation.
  enum class KeyIndex : size_t {
    Parent_Id,
    Child_Id,
    Edge,
    Parent_RHat,
    Parent_RFocal,
    Parent_PHatSister,
    Child_P,
    Child_PHatLeft,
    Child_PHatRight,
  };
  static const size_t KeyIndexCount = 9;
  using KeyIndexIterator =
      EnumIterator<KeyIndex, KeyIndex::Parent_Id, KeyIndex::Child_PHatRight>;
  using KeyIndexMap = EnumArray<KeyIndex, KeyIndexCount, size_t>;
  using KeyIndexMapPair = std::pair<KeyIndexMap, KeyIndexMap>;
  using KeyIndexPairArray = std::array<std::pair<KeyIndex, KeyIndex>, 3>;

  // Translate NNIClade type to corresponding PHat PLV from KeyIndex type.
  using NNIClade = NNIOperation::NNIClade;
  static KeyIndex NNICladeToPHatPLV(NNIClade clade_type);

  // Builds an array containing a mapping of plvs from Pre-NNI to Post-NNI, according
  // to remapping of the NNI clades (sister, right_child, left_child).
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
      const KeyIndexMap &pre_key_idx);

  // Fetches Pre-NNI data to prep Post-NNI for likelihood computation. Method stores
  // intermediate values in the GPEngine temp space (expects GPEngine has already been
  // resized).
  KeyIndexMapPair PassDataFromPreNNIToPostNNIViaReference(
      const NNIOperation &pre_nni, const NNIOperation &post_nni,
      const size_t nni_count = 0, const bool use_unique_temps = false);
  // Fetches Data from Pre-NNI for Post-NNI. Can be used for initial values when
  // moving accepted NNIs from Graft to Host DAG.
  KeyIndexMapPair PassDataFromPreNNIToPostNNIViaCopy(const NNIOperation &pre_nni,
                                                     const NNIOperation &post_nni);

  // Build GPOperation vector for computing NNI Likelihood. For NNIs that are not in
  // the DAG, must create KeyIndexMap through PassDataFromPreNNIToPostNNIViaReference.
  GPOperationVector BuildGPOperationsForNNILikelihood(
      const NNIOperation &nni, const KeyIndexMap &nni_key_idx) const;
  // Build GPOperation vector for computing NNI Likelihood for all adjacent NNIs.
  GPOperationVector BuildGPOperationsForAdjacentNNILikelihoods(
      const bool via_reference = true);
  // Grow engine to handle computing NNI Likelihoods for all adjacent NNIs.
  // Option to grow engine for computing likelihoods via reference or via copy. If
  // computing via reference, option whether to use unique temporaries (for testing
  // and computing in parallel).
  void GrowEngineForAdjacentNNILikelihoods(const bool via_reference = true,
                                           const bool use_unique_temps = false);
  // Performs entire likelihood computation for all Adjacent NNIs.
  // Allocates necessary extra space on GPEngine, generates and processes
  // GPOperationVector, then retrieves results and stores in Scored NNIs.
  void ScoreAdjacentNNIsByLikelihood();

  // ** Set Filtering Scheme

  // Set custom functions for filtering NNIs.
  void SetFilteringCustom(StaticFilterInitFunction init_fn,
                          StaticFilterUpdateFunction pre_update_fn,
                          StaticFilterProcessFunction process_fn,
                          StaticFilterUpdateFunction post_update_fn);
  // Set no filters for NNIs. Accept parameter determines whether to accept or reject
  // all NNIs.
  void SetFilteringNone();
  // Set filter to accept only Post-NNIs which are an improvement over Pre-NNIs.
  void SetFilteringHillclimb();
  // Set filter to compare Pre-NNI likelihood to Post-NNI likelihood.  Accept Post-NNIs
  // whose loss in likelihood relative to its Pre-NNI is less than the loss_tolerance.
  // A loss_tolerance of zero is equivalent to hill-climbing.
  void SetFilteringLossThreshold(const double loss_tolerance);
  // Set filter to have a hard cutoff for accepted Post-NNI likelihood.
  void SetFilteringCutoffThreshold(const double cutoff_threshold);

  // ** DAG & GPEngine Maintenance

  // Initial GPEngine for use with GraftDAG.
  void GPEngineInit();
  // Resize GPEngine for use with GraftDAG, Populate PLVs and Compute Likelihoods for
  // HostDAG.
  void GPEngineUpdate();

  // Add all Accepted NNIs to Main DAG.
  void AddAcceptedNNIsToDAG();
  // Add all Adjacent NNIs to Graft DAG.
  void GraftAdjacentNNIsToDAG();
  // Remove all NNIs from Graft DAG.
  void RemoveAllGraftedNNIsFromDAG();

  // ** NNI Maintenance
  // These maintain NNIs to stay consistent with the state of associated GraftDAG.

  // Freshly synchonizes NNISet to match the current state of its DAG. Wipes old NNI
  // data and finds all all parent/child pairs adjacent to DAG by iterating over all
  // internal edges in the DAG. (For each internal edges, two NNIs are possible.)
  void SyncAdjacentNNIsWithDAG();
  // Updates NNI Set after given parent/child node pair have been added to the DAG.
  // Removes pair from NNI Set and adds adjacent pairs coming from newly created
  // edges.
  void UpdateAdjacentNNIsAfterDAGAddNodePair(const NNIOperation &nni);
  void UpdateAdjacentNNIsAfterDAGAddNodePair(const Bitset &parent_bitset,
                                             const Bitset &child_bitset);
  // Adds all NNIs from all (node_id, other_id) pairs, where other_id's are elements
  // of the adjacent_node_ids vector. is_edge_leafward tells whether node_id is the
  // child or parent. is_edge_edgeclade determines which side of parent the child
  // descends from.
  void AddAllNNIsFromNodeVectorToAdjacentNNIs(const size_t &node_id,
                                              const SizeVector &adjacent_node_ids,
                                              const bool is_edge_on_left,
                                              const bool is_edge_leafward);
  // Based on given input NNIOperation, produces the two possible output NNIOperations
  // and adds those results to the NNI Set (if results are not a member of the DAG or
  // NNI Set).
  void SafeAddOutputNNIsToAdjacentNNIs(const Bitset &parent_bitset,
                                       const Bitset &child_bitset,
                                       const bool is_edge_on_left);

  // Add score to given NNI.
  void AddScoreForNNI(const NNIOperation &nni, const double score);
  // Get score for given NNI.
  double GetScoreForNNI(const NNIOperation &nni) const;

  // Update adjacent NNIs at end of current sweep. If not re-evaluating rejected NNIs,
  // then current adjacent NNIs are removed. Then NNIs that are adjacent to last
  // sweep's accepted NNIs are added.
  void UpdateAdjacentNNIs();
  // Remove all accepted NNIs and optionally save to past NNIs.
  void UpdateAcceptedNNIs();
  // Remove all rejected NNIs and optionally save to past NNIs.
  void UpdateRejectedNNIs();
  // Remove all scored NNIs and optionally save to past NNIs.
  void UpdateScoredNNIs();
  // Reset all NNIs, current and past.
  void ResetAllNNIs();

  // ** Miscellaneous

  NNIOperation FindNNINeighborInDAG(const NNIOperation &nni);

 private:
  // Un-owned reference DAG.
  GPDAG &dag_;
  // For adding temporary NNIs to DAG.
  std::unique_ptr<GraftDAG> graft_dag_;
  // Un-owned reference GPEngine.
  GPEngine &gp_engine_;
  // Tracks node modifications to the DAG.
  Reindexer node_reindexer_;
  // Tracks edge modifications to the DAG.
  Reindexer edge_reindexer_;

  // Set of NNIs to be evaluated, which are a single NNI.
  NNISet adjacent_nnis_;
  // NNIs which have passed the filtering threshold, to be added to the DAG.
  NNISet accepted_nnis_;
  // NNIs which have been accepted in a previous sweep of the search.
  NNISet accepted_past_nnis_;
  // Determines whether past accepted NNIs are tracked.
  bool keep_accepted_past_nnis_ = false;
  // NNIs which have failed the filtering threshold, not to be added to the DAG.
  NNISet rejected_nnis_;
  // NNIs which have been rejected in a previous sweep of the search.
  NNISet rejected_past_nnis_;
  // Determines whether past rejected NNIs are tracked.
  bool keep_rejected_past_nnis_ = false;
  // Determines whether to keep rejected NNIs in Adjacent NNI list for re-evaluation.
  bool reevaluate_rejected_nnis_ = false;

  // Map of adjacent NNIs to their score.
  NNIDoubleMap scored_nnis_;
  // Map of previous rejected NNIs to their score.
  NNIDoubleMap scored_past_nnis_;
  // Determines whether past scored NNIs are tracked.
  bool keep_scored_past_nnis_ = false;

  // DAG Likelihoods.
  std::unique_ptr<EigenVectorXd> dag_likelihoods_ = nullptr;

  // Determines whether branch lengths are optimized when updating DAG likelihoods.
  bool do_optimitize_branch_lengths_ = false;

  // Initial step before running any sweep iterations.
  StaticFilterInitFunction filter_init_fn_ = NoInit;
  // Update step for any filter computation work to be done at start of each sweep.
  StaticFilterUpdateFunction filter_pre_update_fn_ = NoUpdate;
  // Update step for any filter computation work to be done at end of each sweep.
  StaticFilterUpdateFunction filter_post_update_fn_ = NoUpdate;
  // Processing an adjacent NNI to be accepted or rejected.
  StaticFilterProcessFunction filter_process_fn_ = NoFilter<true>;

  // Count number of loops executed by engine.
  size_t sweep_count_ = 0;
};

template <typename DAGType>
NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForNNI(const NNIOperation &nni,
                                                         const DAGType &dag,
                                                         const size_t node_count) {
  // Find NNI nodes.
  const auto parent_id = dag.GetDAGNodeId(nni.GetParent());
  const auto child_id = dag.GetDAGNodeId(nni.GetChild());
  const bool is_left_clade_sister = nni.WhichParentCladeIsFocalClade();
  // Find key indices for NNI.
  KeyIndexMap key_idx_map(NoId);
  key_idx_map[KeyIndex::Parent_Id] = parent_id;
  key_idx_map[KeyIndex::Child_Id] = child_id;
  key_idx_map[KeyIndex::Edge] = dag.GetEdgeIdx(parent_id, child_id);
  key_idx_map[KeyIndex::Parent_RHat] =
      PLVHandler::GetPLVIndex(PLVType::RHat, parent_id, node_count);
  key_idx_map[KeyIndex::Parent_RFocal] = PLVHandler::GetPLVIndex(
      PLVHandler::RPLVType(!is_left_clade_sister), parent_id, node_count);
  key_idx_map[KeyIndex::Parent_PHatSister] = PLVHandler::GetPLVIndex(
      PLVHandler::PPLVType(is_left_clade_sister), parent_id, node_count);
  key_idx_map[KeyIndex::Child_P] =
      PLVHandler::GetPLVIndex(PLVType::P, child_id, node_count);
  key_idx_map[KeyIndex::Child_PHatLeft] =
      PLVHandler::GetPLVIndex(PLVType::PHatLeft, child_id, node_count);
  key_idx_map[KeyIndex::Child_PHatRight] =
      PLVHandler::GetPLVIndex(PLVType::PHatRight, child_id, node_count);

  return key_idx_map;
}

template <typename DAGType>
NNIEngine::KeyIndexMap NNIEngine::BuildKeyIndexMapForPostNNIViaReferencePreNNI(
    const NNIOperation &pre_nni, const NNIOperation &post_nni,
    const NNIEngine::KeyIndexMap &pre_key_idx, const DAGType &dag) {
  // Unpopulated key indices will left as NoId.
  KeyIndexMap post_key_idx(NoId);
  post_key_idx[KeyIndex::Parent_Id] = dag.GetDAGNodeId(post_nni.GetParent());
  post_key_idx[KeyIndex::Child_Id] = dag.GetDAGNodeId(post_nni.GetChild());
  post_key_idx[KeyIndex::Edge] = dag.GetEdgeIdx(post_key_idx[KeyIndex::Parent_Id],
                                                post_key_idx[KeyIndex::Child_Id]);

  // Array for mapping from pre-NNI plvs to post-NNI plvs.
  const auto key_map = BuildKeyIndexTypePairsFromPreNNIToPostNNI(pre_nni, post_nni);

  // Set NNI plvs to their corresponding Pre-NNI plvs.
  post_key_idx[KeyIndex::Parent_RHat] = pre_key_idx[KeyIndex::Parent_RHat];
  for (const auto &[pre_key_type, post_key_type] : key_map) {
    post_key_idx[post_key_type] = pre_key_idx[pre_key_type];
  }

  return post_key_idx;
}
