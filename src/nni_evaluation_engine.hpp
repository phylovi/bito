// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The NNI Evaluation Engine is an interface between the NNI Engine, that performs the
// NNI systematic search over the DAG, finding

#pragma once

#include "gp_dag.hpp"
#include "nni_engine.hpp"

using NNIDoubleMap = std::map<NNIOperation, double>;

class NNIEvaluationEngine {
 public:
  using KeyIndex = NNIEngine::KeyIndex;
  using KeyIndexMap = NNIEngine::KeyIndexMap;
  using KeyIndexMapPair = NNIEngine::KeyIndexMapPair;

  explicit NNIEvaluationEngine(NNIEngine &nni_engine)
      : nni_engine_(&nni_engine),
        dag_(&nni_engine.GetGPDAG()),
        graft_dag_(&nni_engine.GetGraftDAG()) {}

  // Initialize Evaluation Engine.
  virtual void Init();
  // Prepare Evaluation Engine for NNI Engine loop.
  virtual void Prep();
  // Resize Engine for modified DAG.
  virtual void GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                                std::optional<Reindexer> edge_reindexer);
  // Grow engine to handle computing NNIs for all adjacent NNIs.
  // Option to grow engine for computing via reference or via copy. If computing via
  // reference, option whether to use unique temporaries (for testing and computing in
  // parallel).
  virtual void GrowEngineForAdjacentNNILikelihoods(const NNISet &adjacent_nnis,
                                                   const bool via_reference,
                                                   const bool use_unique_temps);
  // Grow Evaluation Engine to account for
  // Compute scores for all NNIs adjacent to current DAG.
  virtual void ScoreAdjacentNNIs(const NNISet &adjacent_nnis);

  // ** Access

  // Get reference NNIEngine.
  NNIEngine &GetNNIEngine() { return *nni_engine_; }
  const NNIEngine &GetNNIEngine() const { return *nni_engine_; }
  // Get reference DAG.
  const GPDAG &GetDAG() const { return *dag_; }
  // Get reference GraftDAG.
  const GraftDAG &GetGraftDAG() const { return *graft_dag_; }
  // Get all Scored NNIs.
  NNIDoubleMap &GetScoredNNIs() { return scored_nnis_; }
  const NNIDoubleMap &GetScoredNNIs() const { return scored_nnis_; }

  // Retrieve Score for given NNI.
  double GetScoreByNNI(const NNIOperation &nni) const {
    return GetScoredNNIs().find(nni)->second;
  }
  // Retrieve Score for given edge in DAG.
  double GetScoreByEdge(const EdgeId edge_id) const {
    auto nni = GetDAG().GetNNI(edge_id);
    return GetScoreByNNI(nni);
  }

 private:
  // Un-owned reference to NNIEngine.
  NNIEngine *nni_engine_ = nullptr;
  // Un-owned reference to DAG.
  const GPDAG *dag_ = nullptr;
  // Un-owned reference to GraftDAG.
  const GraftDAG *graft_dag_ = nullptr;
  // Scored NNIs.
  NNIDoubleMap scored_nnis_;
};

class NNIEvaluationEngineViaGP : public NNIEvaluationEngine {
 public:
  NNIEvaluationEngineViaGP(NNIEngine &nni_engine, GPEngine &gp_engine)
      : NNIEvaluationEngine(nni_engine), gp_engine_(&gp_engine) {}

  // Initialize Evaluation Engine.
  virtual void Init();
  // Prepare Engine for NNI Engine loop.
  virtual void Prep();
  // Resize Engine for modified DAG.
  virtual void GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                                std::optional<Reindexer> edge_reindexer);
  // Fetches Pre-NNI data to prep Post-NNI for likelihood computation. Method stores
  // intermediate values in the GPEngine temp space (expects GPEngine has already been
  // resized).
  virtual void GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                         const bool via_reference,
                                         const bool use_unique_temps);
  // Compute scores for all NNIs adjacent to current DAG.
  virtual void ScoreAdjacentNNIs(const NNISet &adjacent_nnis);

  // ** Helpers

  // Fetches Data from Pre-NNI for Post-NNI. Can be used for initial values when moving
  // accepted NNIs from Graft to Host DAG.
  GPOperationVector BuildGPOperationsForAdjacentNNILikelihoods(
      const NNISet &adjacent_nnis, const bool via_reference);
  // Build GPOperation vector for computing NNI Likelihood. For NNIs that are not in
  // the DAG, must create KeyIndexMap through
  // PassGPEngineDataFromPreNNIToPostNNIViaReference.
  GPOperationVector BuildGPOperationsForNNILikelihood(
      const NNIOperation &nni, const NNIEngine::KeyIndexMap &nni_key_idx) const;
  // Build GPOperation vector for computing NNI Likelihood for all adjacent NNIs.
  GPOperationVector BuildGPOperationsForAdjacentNNILikelihoods(
      const NNISet &adjacent_nnis, const bool via_reference);
  // Grow engine to handle computing NNI Likelihoods for all adjacent NNIs.
  // Option to grow engine for computing likelihoods via reference or via copy. If
  // computing via reference, option whether to use unique temporaries (for testing and
  // computing in parallel).
  void GrowGPEngineForAdjacentNNILikelihoods(const NNISet &adjacent_nnis,
                                             const bool via_reference,
                                             const bool use_unique_temps);
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

  // ** Access

  GPEngine &GetGPEngine() { return *gp_engine_; }

 private:
  GPEngine *gp_engine_ = nullptr;
};

class NNIEvaluationEngineViaTP : public NNIEvaluationEngine {
 public:
  using KeyIndex = NNIEngine::KeyIndex;
  using KeyIndexMap = NNIEngine::KeyIndexMap;

  NNIEvaluationEngineViaTP(NNIEngine &nni_engine, TPEngine &tp_engine)
      : NNIEvaluationEngine(nni_engine), tp_engine_(&tp_engine) {}

  // Initialize Evaluation Engine.
  virtual void Init();
  // Prepare Engine for NNI Engine loop.
  virtual void Prep();
  // Resize Engine for modified DAG.
  virtual void GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                                std::optional<Reindexer> edge_reindexer);
  // Fetches Pre-NNI data to prep Post-NNI for likelihood computation. Method stores
  // intermediate values in the GPEngine temp space (expects GPEngine has already been
  // resized).
  virtual void GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                         const bool via_reference,
                                         const bool use_unique_temps);
  // Compute scores for all NNIs adjacent to current DAG.
  virtual void ScoreAdjacentNNIs(const NNISet &adjacent_nnis);

  // ** Access

  TPEngine &GetTPEngine() { return *tp_engine_; }

 private:
  TPEngine *tp_engine_;
};
