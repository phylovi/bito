// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// NNI Computation Engine
// This is an interface for computation engines for scoring edges within the DAG and
// adjacent to the DAG.

class NNIComputationEngine {
 public:
  // Performs entire TP likelihood computation for all Adjacent NNIs.
  // Allocates necessary extra space on engine, computes and retrieves results and
  // stores in Scored NNIs.
  virtual void ScoreAdjacentNNIs();
  // Initial engine for use with GraftDAG.
  virtual void InitEngine();
  // Populate PLVs for quick lookup of scores.
  virtual void PrepEngine();
};
