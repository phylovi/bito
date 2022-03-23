// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "subsplit_dag.hpp"
#include "node.hpp"
#include "mersenne_twister.hpp"

class TopologySampler {
 public:

  Node::NodePtr Sample(SubsplitDAGNode node,
    SubsplitDAG& dag,
    EigenConstVectorXdRef normalized_sbn_parameters,
    EigenConstVectorXdRef inverted_probabilities);

  void SetSeed(uint64_t seed);

 private:

  struct SamplingSession {
    SubsplitDAG& dag_;
    EigenConstVectorXdRef normalized_sbn_parameters_;
    EigenConstVectorXdRef inverted_probabilities_;
    SubsplitDAGStorage result_;
  };

  void VisitNode(SamplingSession& session, SubsplitDAGNode node, Direction direction, Clade clade);
  void SampleRootward(SamplingSession& session, SubsplitDAGNode node);
  void SampleLeafward(SamplingSession& session, SubsplitDAGNode node, Clade clade);
  std::pair<SubsplitDAGNode, ConstLineView> SampleParent(SamplingSession& session, ConstNeighborsView left, ConstNeighborsView right);
  std::pair<SubsplitDAGNode, ConstLineView> SampleChild(SamplingSession& session, ConstNeighborsView neighbors);
  Node::NodePtr BuildTree(SamplingSession& session, const DAGVertex& node);

  MersenneTwister mersenne_twister_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("TopologySampler") {
  Driver driver;
  auto tree_collection = RootedTreeCollection::OfTreeCollection(
      driver.ParseNewickFile("data/five_taxon_rooted.nwk"));
  SubsplitDAG dag(tree_collection);

  MersenneTwister mersenne_twister;
  std::uniform_int_distribution<> distribution(0, dag.TaxonCount());
  auto leaf = dag.GetDAGNode(distribution(mersenne_twister.GetGenerator()));

  EigenVectorXd normalized_sbn_parameters = dag.BuildUniformOnTopologicalSupportPrior();
  EigenVectorXd node_probabilities =
      dag.UnconditionalNodeProbabilities(normalized_sbn_parameters);
  EigenVectorXd inverted_probabilities =
      dag.InvertedGPCSPProbabilities(normalized_sbn_parameters, node_probabilities);

  TopologySampler sampler;
  auto tree = sampler.Sample(leaf, dag, normalized_sbn_parameters, inverted_probabilities);
  std::ignore = tree;
}

#endif  // DOCTEST_LIBRARY_INCLUDED
