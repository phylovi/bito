#pragma once

#include "../src/topology_sampler.hpp"

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

  SubsplitDAGNode origin = dag.GetDAGNode(NodeId(5));

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

  SubsplitDAGNode origin = dag.GetDAGNode(NodeId(5));

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

  for (auto& [tree, count] : counts) {
    const double observed = static_cast<double>(count) / iterations;
    CHECK_LT(fabs(observed - expected[tree]), 5e-2);
  }
}

#endif  // DOCTEST_LIBRARY_INCLUDED
