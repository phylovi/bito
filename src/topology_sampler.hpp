// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "subsplit_dag.hpp"
#include "node.hpp"
#include "mersenne_twister.hpp"

class TopologySampler {
 public:
  Node::NodePtr Sample(SubsplitDAGNode node, SubsplitDAG& dag,
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

  void VisitNode(SamplingSession& session, SubsplitDAGNode node, Direction direction,
                 Clade clade);
  void SampleRootward(SamplingSession& session, SubsplitDAGNode node);
  void SampleLeafward(SamplingSession& session, SubsplitDAGNode node, Clade clade);
  std::pair<SubsplitDAGNode, ConstLineView> SampleParent(SamplingSession& session,
                                                         ConstNeighborsView left,
                                                         ConstNeighborsView right);
  std::pair<SubsplitDAGNode, ConstLineView> SampleChild(SamplingSession& session,
                                                        ConstNeighborsView neighbors);
  Node::NodePtr BuildTree(SamplingSession& session, const DAGVertex& node);

  MersenneTwister mersenne_twister_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("TopologySampler") {
  Driver driver;
  auto tree_collection = RootedTreeCollection::OfTreeCollection(
      driver.ParseNewickFile("data/five_taxon_rooted_more_2.nwk"));
  SubsplitDAG dag(tree_collection);

  EigenVectorXd normalized_sbn_parameters = dag.BuildUniformOnTopologicalSupportPrior();
  EigenVectorXd node_probabilities =
      dag.UnconditionalNodeProbabilities(normalized_sbn_parameters);
  EigenVectorXd inverted_probabilities =
      dag.InvertedGPCSPProbabilities(normalized_sbn_parameters, node_probabilities);

  SubsplitDAGNode origin = dag.GetDAGNode(5);

  TopologySampler sampler;
  std::map<std::string, size_t> counts;
  const size_t iterations = 10000;

  for (size_t i = 0; i < iterations; ++i) {
    auto tree =
        sampler.Sample(origin, dag, normalized_sbn_parameters, inverted_probabilities);
    ++counts[tree->Newick([](const Node* node) {
      if (!node->IsLeaf()) return std::string();
      return std::string("x") + std::to_string(node->Id());
    })];
  }

  for (auto& i : counts) {
    const double observed = static_cast<double>(i.second) / iterations;
    const double expected = 1.0 / 3.0;
    std::cout << i.first << " : " << observed << "\n";
    CHECK_LT(fabs(observed - expected), 5e-2);
  }
}

TEST_CASE("TopologySampler: Non-uniform prior") {
  Driver driver;
  auto tree_collection = RootedTreeCollection::OfTreeCollection(
      driver.ParseNewickFile("data/five_taxon_rooted_more_2.nwk"));
  SubsplitDAG dag(tree_collection);

  std::vector<double> params{0.5, 0.3, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0,
                             0.8, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  EigenVectorXd normalized_sbn_parameters =
      EigenVectorXd::Map(params.data(), params.size());
  EigenVectorXd node_probabilities =
      dag.UnconditionalNodeProbabilities(normalized_sbn_parameters);
  EigenVectorXd inverted_probabilities =
      dag.InvertedGPCSPProbabilities(normalized_sbn_parameters, node_probabilities);

  SubsplitDAGNode origin = dag.GetDAGNode(5);

  TopologySampler sampler;
  std::map<std::string, size_t> counts;
  std::map<std::string, double> expected = {
      {"((((x0,x1),x2),(x3,x4)));", 0.312},
      {"(((x0,x1),(x2,(x3,x4))));", 0.52},
      {"((x0,(x1,(x2,(x3,x4)))));", 0.1666},
  };
  const size_t iterations = 10000;

  for (size_t i = 0; i < iterations; ++i) {
    auto tree =
        sampler.Sample(origin, dag, normalized_sbn_parameters, inverted_probabilities);
    ++counts[tree->Newick([](const Node* node) {
      if (!node->IsLeaf()) return std::string();
      return std::string("x") + std::to_string(node->Id());
    })];
  }

  for (auto& i : counts) {
    const double observed = static_cast<double>(i.second) / iterations;
    CHECK_LT(fabs(observed - expected[i.first]), 5e-2);
  }
}

#endif  // DOCTEST_LIBRARY_INCLUDED
