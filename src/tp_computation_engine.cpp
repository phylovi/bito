// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// TP Computation Engine

#pragma once

#include "tp_engine.hpp"

class TPComputationEngine {
 public:
  // Initialize Computation Engine.
  virtual void Init();
  // Update the Computation Engine after adding Node Pairs to the DAG.
  virtual void UpdateAfterDAGNodePair(const NNIOperation &post_nni,
                                      const NNIOperation &pre_nni,
                                      std::optional<size_t> new_tree_id);
  // Get the Top Tree from the DAG with the given edge.
  virtual double GetTopTreeLikelihoodWithEdge(const EdgeId edge_id);
};

class TPComputationEngineViaLikelihood : public TPComputationEngine {
 public:
  // Initialize Computation Engine.
  virtual void Init();
  // Update the Computation Engine after adding Node Pairs to the DAG.
  virtual void UpdateAfterDAGNodePair(const NNIOperation &post_nni,
                                      const NNIOperation &pre_nni,
                                      std::optional<size_t> new_tree_id);
  // Get the Top Tree from the DAG with the given edge.
  virtual double GetTopTreeLikelihoodWithEdge(const EdgeId edge_id);
};

class TPComputationEngineViaParsimony : public TPComputationEngine {
 public:
  // Initialize Computation Engine.
  virtual void Init();
  // Update the Computation Engine after adding Node Pairs to the DAG.
  virtual void UpdateAfterDAGNodePair(const NNIOperation &post_nni,
                                      const NNIOperation &pre_nni,
                                      std::optional<size_t> new_tree_id);
  // Get the Top Tree from the DAG with the given edge.
  virtual double GetTopTreeLikelihoodWithEdge(const EdgeId edge_id);
};
