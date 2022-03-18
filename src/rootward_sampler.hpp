// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "subsplit_dag.hpp"
#include "node.hpp"
#include "mersenne_twister.hpp"

class RootwardSampler {
 public:
  Node::NodePtr SampleTopology(SubsplitDAGNode node,
                               EigenConstVectorXdRef inverted_probabilities);

  void SetSeed(uint64_t seed);

 private:
  void SampleTopology(SubsplitDAGNode node, std::vector<size_t>& edges,
                      EigenConstVectorXdRef inverted_probabilities);

  SubsplitDAGNode SampleParent(ConstNeighborsView left, ConstNeighborsView right,
                               EigenConstVectorXdRef inverted_probabilities);

  MersenneTwister mersenne_twister_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("RootwardSampler") {
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

  RootwardSampler sampler;
  auto tree = sampler.SampleTopology(leaf, inverted_probabilities);
  std::ignore = tree;
}

#endif  // DOCTEST_LIBRARY_INCLUDED
