// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// TP Computation Engine

#pragma once

#include "tp_engine.hpp"

class TPComputationEngine {
  //
  virtual void Initialize();
  //
  virtual void UpdateAfterDAGNodePair(const NNIOperation &post_nni,
                                      const NNIOperation &pre_nni,
                                      std::optional<size_t> new_tree_id);
  //
  virtual double GetTopTreeLikelihoodWithEdge(const EdgeId edge_id);
}
