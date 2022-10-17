// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The NNI Evaluation Engine is an interface between the NNI Engine, that performs the
// NNI systematic search over the DAG, finding

#pragma once

#include "gp_dag.hpp"
#include "nni_engine.hpp"

class NNIEvaluationEngine {
 public:
  explicit NNIEvaluationEngine(GPDAG &dag) : dag_(&dag) {}

  // Initialize Evaluation Engine.
  virtual void Init();
  // Prepare Evaluation Engine for NNI Engine loop.
  virtual void Prep();
  // Compute scores for all NNIs adjacent to current DAG.
  virtual void ScoreAdjacentNNIs();

  // ** Access

  // Get reference DAG.
  GPDAG &GetDAG() { return *dag_; }
  // Get all Scored NNIs.
  NNIDoubleMap &GetScoredNNIs() { return *scored_nnis_; }

  // Retrieve Score for given NNI.
  double GetScoreByNNI(const NNIOperation &nni) const {
    return GetScoredNNIs().find(nni)->second;
  }
  // Retrieve Score for given edge in DAG.
  double GetScoreByEdge(const EdgeId edge_id) const {}

  GPDAG *dag_;
  NNIDoubleMap *scored_nnis_;
};

class NNIEvaluationEngineViaGP : public NNIEvaluationEngine {
 public:
  explicit NNIEvaluationEngineViaGP(GPDAG &dag, GPEngine &gp_engine)
      : NNIEvaluationEngine(dag), gp_engine_(&gp_engine) {}

  // Initialize Evaluation Engine.
  virtual void Init() {
    GetGPEngine().GrowPLVs(GetDAG().NodeCountWithoutDAGRoot());
    GetGPEngine().GrowGPCSPs(GetDAG().EdgeCountWithLeafSubsplits());
    GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());
    GetGPEngine().ProcessOperations(GetDAG().ComputeLikelihoods());
  }
  // Prepare Evaluation Engine for NNI Engine loop.
  virtual void Prep() {
    GetGPEngine().ProcessOperations(GetDAG().PopulatePLVs());
    GetGPEngine().ProcessOperations(GetDAG().ComputeLikelihoods());
  }
  // Compute scores for all NNIs adjacent to current DAG.
  virtual void ScoreAdjacentNNIs() {}

  // ** Access

  GPEngine &GetGPEngine() { return *gp_engine_; }

  GPEngine *gp_engine_ = nullptr;
};

class NNIEvaluationEngineViaTP : public NNIEvaluationEngine {
 public:
  NNIEvaluationEngineViaGP(GPDAG &dag, GPEngine &tp_engine)
      : NNIEvaluationEngine(dag), tp_engine_(&tp_engine) {}

  // Initialize Evaluation Engine.
  virtual void Init() {
    GetTPEngine().GrowNodeData(GetDAG().NodeCount());
    GetTPEngine().GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits());
  }
  // Prepare Evaluation Engine for NNI Engine loop.
  virtual void Prep() {}
  // Compute scores for all NNIs adjacent to current DAG.
  virtual void ScoreAdjacentNNIs() {}

  TPEngine *tp_engine_;
};
