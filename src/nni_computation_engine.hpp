// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// NNI Computation Engine

#pragma once

#include "nni_engine.hpp"

class NNIComputationEngine {
  //
  virtual void Init();
  // Prepare for
  virtual void Prep();
  //
  virtual void ScoreAdjacentNNIs();
  //
  virtual double GetScoreByNNI(const NNIOperation &nni) const;
  //
  virtual double GetScoreByEdge(const EdgeId edge_id) const;
};
