// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

// A class that samples one tree topology per call to the Sample() function.
// The parameters are as follows:
//   node - a node to start sampling for
//   dag - the SubsplitDAG that owns "node"
//   normalized_sbn_parameters - edge probabilities for leafward sampling
//   inverted_probabilities - edge probabilities for rootward sampling

#pragma once

#include "subsplit_dag.hpp"
#include "node.hpp"
#include "mersenne_twister.hpp"

class TopologySampler {
 public:
  // Sample a single tree from the DAG according to the provided edge
  // probabilities.
  Node::NodePtr Sample(SubsplitDAGNode node, SubsplitDAG& dag,
                       EigenConstVectorXdRef normalized_sbn_parameters,
                       EigenConstVectorXdRef inverted_probabilities);

  // Set a seed value for the internal random number generator.
  void SetSeed(uint64_t seed);

 private:
  struct SamplingSession {
    SubsplitDAG& dag_;
    EigenConstVectorXdRef normalized_sbn_parameters_;
    EigenConstVectorXdRef inverted_probabilities_;
    SubsplitDAGStorage result_;
  };

  // Called for each newly sampled node. Direction and clade are pointing to the
  // previously visited node.
  void VisitNode(SamplingSession& session, SubsplitDAGNode node, Direction direction,
                 SubsplitClade clade);
  // Continue sampling in the rootward direction from `node`.
  void SampleRootward(SamplingSession& session, SubsplitDAGNode node);
  // Continue sampling in the leafward direction and specified `clade` from `node`.
  void SampleLeafward(SamplingSession& session, SubsplitDAGNode node,
                      SubsplitClade clade);
  // Choose a parent node (and return it with the corresponding edge) according to
  // the values in `inverted_probabilities`.
  std::pair<SubsplitDAGNode, ConstLineView> SampleParentNodeAndEdge(
      SamplingSession& session, ConstNeighborsView left, ConstNeighborsView right);
  // Choose a child node among the `neighbors` clade according to the values
  // in `normalized_sbn_parameters`.
  std::pair<SubsplitDAGNode, ConstLineView> SampleChildNodeAndEdge(
      SamplingSession& session, ConstNeighborsView neighbors);
  // Construct a Node topology from a successful sampling. Recursion is started
  // by passing the root node.
  Node::NodePtr BuildTree(SamplingSession& session, const DAGVertex& node);

  MersenneTwister mersenne_twister_;
};
