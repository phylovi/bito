// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The NNI Evaluation Engine is an interfaces for the NNI Engine, that facilitates
// different methods for scoring NNIs. The NNIEngine procedurally finds all NNIs
// adjacent to the current DAG, then passes these NNIs to the helper NNIEvalEngine
// class for scoring. Currently supports scoring by Generalized Pruning and Top Pruning
// (using Likelihoods or Parsimonies).

#pragma once

#include "gp_dag.hpp"
#include "graft_dag.hpp"
#include "nni_engine_key_index.hpp"
#include "dag_branch_handler.hpp"
#include "gp_engine.hpp"
#include "tp_engine.hpp"

using NNIDoubleMap = std::map<NNIOperation, double>;

// Forward declaration.
class NNIEngine;

// NNIEngine helper for evaluation NNIs -- base class.
class NNIEvalEngine {
 public:
  using KeyIndex = NNIEngineKeyIndex;
  using KeyIndexPairArray = NNIEngineKeyIndexPairArray;
  using KeyIndexMap = NNIEngineKeyIndexMap;
  using KeyIndexMapPair = NNIEngineKeyIndexMapPair;

  explicit NNIEvalEngine(NNIEngine &nni_engine);

  virtual ~NNIEvalEngine() {}

  // ** Maintenance

  // Initialize Evaluation Engine.
  virtual void Init() { Failwith("Pure virtual function call."); }
  // Prepare Evaluation Engine for NNI Engine loop.
  virtual void Prep() { Failwith("Pure virtual function call."); }
  // Resize Engine for modified DAG.
  virtual void GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                                std::optional<Reindexer> edge_reindexer) {
    Failwith("Pure virtual function call.");
  }
  // Grow engine to handle computing NNIs for all adjacent NNIs.
  // Option to grow engine for computing via reference or via copy. If computing via
  // reference, option whether to use unique temporaries (for testing and computing in
  // parallel).
  virtual void GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                         const bool via_reference = true,
                                         const bool use_unique_temps = true) {
    Failwith("Pure virtual function call.");
  }
  // Update engine after modifying DAG (adding nodes and edges).
  virtual void UpdateEngineAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &pre_nni_to_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer) {
    Failwith("Pure virtual function call.");
  }

  // ** Scoring

  // Grow Evaluation Engine to account for
  // Compute scores for all NNIs adjacent to current DAG.
  virtual void ScoreAdjacentNNIs(const NNISet &adjacent_nnis) {
    Failwith("Pure virtual function call.");
  }
  // Score NNI currently in DAG. Expects engine has been prepped and updated after
  // modifying DAG.
  virtual double ScoreInternalNNIByNNI(const NNIOperation &nni) const {
    Failwith("Pure virtual function call.");
  }
  virtual double ScoreInternalNNIByEdge(const EdgeId &edge_id) const {
    Failwith("Pure virtual function call.");
  }
  // Get the number of spare nodes needed per proposed NNI.
  virtual size_t GetSpareNodesPerNNI() const {
    Failwith("Pure virtual function call.");
  }
  // Get the number of spare edges needed per proposed NNI.
  virtual size_t GetSpareEdgesPerNNI() const {
    Failwith("Pure virtual function call.");
  }

  // ** Access

  // Get reference NNIEngine.
  const NNIEngine &GetNNIEngine() const {
    Assert(nni_engine_ != nullptr, "DAG is not set.");
    return *nni_engine_;
  }
  // Get reference DAG.
  const GPDAG &GetDAG() const {
    Assert(dag_ != nullptr, "DAG is not set.");
    return *dag_;
  }
  // Get reference GraftDAG.
  const GraftDAG &GetGraftDAG() const {
    Assert(graft_dag_ != nullptr, "GraftDAG is not set.");
    return *graft_dag_;
  }
  // Get all Scored NNIs.
  const NNIDoubleMap &GetScoredNNIs() const { return scored_nnis_; }
  // Get DAG branch handler.
  const DAGBranchHandler &GetDAGBranchHandler() const {
    Failwith("Pure virtual function call.");
  }

  // Retrieve Score for given NNI.
  double GetScoreByNNI(const NNIOperation &nni) const;
  // Retrieve Score for given edge in DAG.
  double GetScoreByEdge(const EdgeId edge_id) const;
  // Find maximum score in DAG.
  double GetMaxScore() const;
  // Find minimum score in DAG.
  double GetMinScore() const;

  // Determines whether to optimize edges on initialization.
  bool IsOptimizeOnInit() const { return optimize_on_init_; }
  void SetOptimizeOnInit(const bool optimize_on_init) {
    optimize_on_init_ = optimize_on_init;
  }
  // Determines whether new edges are optimized.
  bool IsOptimizeNewEdges() const { return optimize_new_edges_; }
  void SetOptimizeNewEdges(const bool optimize_new_edges) {
    optimize_new_edges_ = optimize_new_edges;
  }
  // Determines whether new edges are initialized by referencing pre-NNI edges.
  bool IsCopyNewEdges() const { return copy_new_edges_; }
  void SetCopyNewEdges(const bool copy_new_edges) { copy_new_edges_ = copy_new_edges; }
  // Determines maximum number of iterations to perform when optimizing edges.
  size_t GetOptimizationMaxIteration() const { return optimize_max_iter_; }
  void SetOptimizationMaxIteration(const size_t optimize_max_iter) {
    optimize_max_iter_ = optimize_max_iter;
  }

 protected:
  // Get reference NNIEngine.
  NNIEngine &GetNNIEngine() {
    Assert(nni_engine_ != nullptr, "NNIEngine is not set.");
    return *nni_engine_;
  }
  // Get reference DAG.
  GPDAG &GetDAG() {
    Assert(dag_ != nullptr, "DAG is not set.");
    return *dag_;
  }
  // Get reference GraftDAG.
  GraftDAG &GetGraftDAG() {
    Assert(graft_dag_ != nullptr, "GraftDAG is not set.");
    return *graft_dag_;
  }
  // Get all Scored NNIs.
  NNIDoubleMap &GetScoredNNIs() { return scored_nnis_; }

  // Un-owned reference to NNIEngine.
  NNIEngine *nni_engine_ = nullptr;
  // Un-owned reference to DAG.
  GPDAG *dag_ = nullptr;
  // Un-owned reference to GraftDAG.
  GraftDAG *graft_dag_ = nullptr;
  // Scored NNIs.
  NNIDoubleMap scored_nnis_;

  // Whether to optimize all edges in DAG on initialization.
  bool optimize_on_init_ = true;
  // Whether new branches are initialized by referencing pre-NNI lengths.
  bool copy_new_edges_ = true;
  // Whether new branches are optimized.
  bool optimize_new_edges_ = true;
  // Number of optimization iterations.
  size_t optimize_max_iter_ = 10;
};

// NNIEngine helper for evaluating NNIs by using Generalized Pruning.  Calls GPEngine
// for functionality.  See gp_engine.hpp for more details.
class NNIEvalEngineViaGP : public NNIEvalEngine {
 public:
  NNIEvalEngineViaGP(NNIEngine &nni_engine, GPEngine &gp_engine);

  // ** Maintenance

  // Initialize EvalEngine.
  void Init() override;
  // Prepare EvalEngine for NNIEngine loop.
  void Prep() override;
  // Resize EvalEngine for modified DAG.
  void GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                        std::optional<Reindexer> edge_reindexer) override;
  // Fetches Pre-NNI data to prep Post-NNI for score computation. Method stores
  // intermediate values in the GPEngine temp space (expects GPEngine has already been
  // resized).
  void GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                 const bool via_reference = true,
                                 const bool use_unique_temps = true) override;
  // Update engine after modifying DAG (adding nodes and edges).
  virtual void UpdateEngineAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &pre_nni_to_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer);

  // ** Scoring

  // Compute scores for all NNIs adjacent to current DAG.
  void ScoreAdjacentNNIs(const NNISet &adjacent_nnis) override;
  // Score NNI currently in DAG. Expects engine has been prepped and updated after
  // modifying DAG.
  double ScoreInternalNNIByNNI(const NNIOperation &nni) const override;
  double ScoreInternalNNIByEdge(const EdgeId &edge_id) const override;
  // Get the number of spare nodes needed per proposed NNI.
  size_t GetSpareNodesPerNNI() const override;
  // Get the number of spare edges needed per proposed NNI.
  size_t GetSpareEdgesPerNNI() const override;

  // ** Helpers

  // Compute GP likelihood for all adjacent NNIs.
  void ComputeAdjacentNNILikelihoods(const NNISet &adjacent_nnis,
                                     const bool via_reference);
  // Build GPOperations for computing GP likelihood for a single adjacent NNI. Returns
  // likelihood and last NNI's offset (number of edges adjacent to NNI).
  std::pair<double, size_t> ComputeAdjacentNNILikelihood(const NNIOperation &nni,
                                                         const size_t offset = 0);

  template <typename T>
  struct AdjEdges {
    std::vector<T> parents;
    std::vector<T> sisters;
    T central;
    std::vector<T> leftchildren;
    std::vector<T> rightchildren;
  };
  using AdjEdgeIds = AdjEdges<EdgeId>;
  using AdjEdgePCSPs = AdjEdges<Bitset>;

  template <typename T>
  struct AdjNodes {
    std::vector<T> grand_parents;
    std::vector<T> grand_sisters;
    T parent_focal;
    T parent_sister;
    T child_left;
    T child_right;
    std::vector<T> grand_leftchildren;
    std::vector<T> grand_rightchildren;
  };
  using AdjNodeIds = AdjNodes<NodeId>;
  using AdjNodeSubsplits = AdjNodes<Bitset>;

  struct AdjPVIds {
    PVId parent_p;
    PVId parent_phatfocal;
    PVId parent_phatsister;
    PVId parent_rhat;
    PVId parent_rfocal;
    PVId parent_rsister;
    PVId child_p;
    PVId child_phatleft;
    PVId child_phatright;
    PVId child_rhat;
    PVId child_rleft;
    PVId child_rright;
  };

  using AdjNodeAndEdgeIds = std::pair<AdjNodeIds, AdjEdgeIds>;
  // Get nodes and edges adjacent to NNI in DAG.
  AdjNodeAndEdgeIds GetAdjNodeAndEdgeIds(const NNIOperation &nni) const;
  // Assign temporary nodes according to pre-NNI in DAG.
  AdjNodeAndEdgeIds GetMappedAdjNodeIdsAndTempAdjEdgeIds(
      const NNIOperation &pre_nni, const NNIOperation &nni,
      const bool copy_branch_lengths = false);
  // Assign temporary nodes according to pre-NNI in DAG.
  AdjEdgeIds GetTempAdjEdgeIds(const AdjNodeIds &node_ids);
  // Assign temporary PVIds for new NNI.
  AdjPVIds GetTempAdjPVIds();

  std::set<EdgeId> BuildSetOfEdgeIdsAdjacentToNNI(const NNIOperation &nni) const;
  std::set<Bitset> BuildSetOfPCSPsAdjacentToNNI(const NNIOperation &nni) const;

  // Copies branch length data from pre-NNI to post-NNI before optimization.
  void CopyGPEngineDataAfterAddingNNI(const NNIOperation &pre_nni,
                                      const NNIOperation &post_nni);

  // Grow engine to handle computing NNI Likelihoods for all adjacent NNIs.
  // Option to grow engine for computing likelihoods via reference or via copy. If
  // computing via reference, option whether to use unique temporaries (for testing and
  // computing in parallel).
  void GrowGPEngineForAdjacentNNILikelihoods(const NNISet &adjacent_nnis,
                                             const bool via_reference = true,
                                             const bool use_unique_temps = true);

  // Fetches Data from Pre-NNI for Post-NNI. Can be used for initial values when moving
  // accepted NNIs from Graft to Host DAG.
  KeyIndexMapPair PassGPEngineDataFromPreNNIToPostNNIViaCopy(
      const NNIOperation &pre_nni, const NNIOperation &post_nni);
  // Fetches Pre-NNI data to prep Post-NNI for likelihood computation. Method stores
  // intermediate values in the GPEngine temp space (expects GPEngine has already been
  // resized).
  KeyIndexMapPair PassGPEngineDataFromPreNNIToPostNNIViaReference(
      const NNIOperation &pre_nni, const NNIOperation &post_nni, const size_t nni_count,
      const bool use_unique_temps);

  // ** Branch Length Optimization

  // Optimize all branch lengths over DAG.
  void BranchLengthOptimization();
  // Optimize only given branch lengths over DAG.
  void BranchLengthOptimization(const std::set<EdgeId> &edges_to_optimize);
  // Optimize only edges adjacent to an NNI.
  void NNIBranchLengthOptimization(const NNIOperation &nni,
                                   const std::set<EdgeId> &new_edge_ids);
  void NNIBranchLengthOptimization(const NNIOperation &nni);

  // ** Access

  // Un-owned reference to GPEngine.
  GPEngine &GetGPEngine() { return *gp_engine_; }
  const GPEngine &GetGPEngine() const { return *gp_engine_; }
  // Get DAG branch handler.
  const DAGBranchHandler &GetDAGBranchHandler() const {
    return GetGPEngine().GetBranchLengthHandler();
  }

 protected:
  // Unowned reference to GPEngine.
  GPEngine *gp_engine_ = nullptr;

  // Spare work space for NNI search.
  size_t spare_nodes_per_nni_ = 2;
  size_t spare_edges_per_nni_ = 5;

  // Whether to use uniform SBN parameters
  bool use_null_priors_ = false;
};

// NNIEngine helper for evaluating NNIs by using Top Pruning.  Calls TPEngine
// for functionality.  See tp_engine.hpp for more details.
class NNIEvalEngineViaTP : public NNIEvalEngine {
 public:
  NNIEvalEngineViaTP(NNIEngine &nni_engine, TPEngine &tp_engine)
      : NNIEvalEngine(nni_engine), tp_engine_(&tp_engine) {}

  // ** Maintenance

  // Initialize Evaluation Engine.
  void Init() override;
  // Prepare Engine for NNI Engine loop.
  void Prep() override;
  // Resize Engine for modified DAG.
  void GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                        std::optional<Reindexer> edge_reindexer) override;
  // Fetches Pre-NNI data to prep Post-NNI for likelihood computation. Method stores
  // intermediate values in the GPEngine temp space (expects GPEngine has already been
  // resized).
  void GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                 const bool via_reference = true,
                                 const bool use_unique_temps = true) override;
  // Update engine after modifying DAG (adding nodes and edges).
  void UpdateEngineAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &pre_nni_to_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer) override;

  // ** Scoring

  // Compute scores for all NNIs adjacent to current DAG.
  void ScoreAdjacentNNIs(const NNISet &adjacent_nnis) override;
  // Score NNI currently in DAG. Expects engine has been prepped and updated after
  // modifying DAG.
  double ScoreInternalNNIByNNI(const NNIOperation &nni) const override;
  double ScoreInternalNNIByEdge(const EdgeId &edge_id) const override;
  // Get the number of spare nodes needed per proposed NNI.
  size_t GetSpareNodesPerNNI() const override;
  // Get the number of spare edges needed per proposed NNI.
  size_t GetSpareEdgesPerNNI() const override;

  // ** Access

  TPEngine &GetTPEngine() { return *tp_engine_; }
  const TPEngine &GetTPEngine() const { return *tp_engine_; }
  // Get DAG branch handler.
  const DAGBranchHandler &GetDAGBranchHandler() const {
    return GetTPEngine().GetDAGBranchHandler();
  }

 protected:
  TPEngine *tp_engine_;
};
