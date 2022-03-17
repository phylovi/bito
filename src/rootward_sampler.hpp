// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "subsplit_dag.hpp"
#include "node.hpp"
#include "mersenne_twister.hpp"

class RootwardSampler {
public:

  Node::NodePtr SampleTopology(SubsplitDAGNode node) {
    std::vector<size_t> edges;
    SampleTopology(node, edges);
    Node::NodePtr tree = Node::Leaf(node.Id());
    for (auto i : edges) {
      tree = Node::Join({tree}, i);
    }
    return tree;
  }

  SubsplitDAGNode RandomLeaf(const SubsplitDAG& dag) {
    std::uniform_int_distribution<> distribution(0, dag.TaxonCount());
    auto i = static_cast<size_t>(distribution(mersenne_twister_.GetGenerator()));
    return dag.GetDAGNode(i);
  }

private:

  void SampleTopology(SubsplitDAGNode node, std::vector<size_t>& edges) {
    auto left = node.GetLeftRootward();
    auto right = node.GetRightRootward();
    if (left.empty() && right.empty()) {
      // reached root
      return;
    }
    auto parent = SampleParent(left, right);
    edges.push_back(parent.Id());
    SampleTopology(parent, edges);
  }

  SubsplitDAGNode SampleParent(ConstNeighborsView left, ConstNeighborsView right) {
    std::uniform_int_distribution<> distribution(0, left.size() + right.size() - 1);
    auto i = static_cast<size_t>(distribution(mersenne_twister_.GetGenerator()));
    if (i < left.size()) {
        auto parent = left.begin();
        std::advance(parent, i);
        return parent.GetNode();
    }
    auto parent = right.begin();
    std::advance(parent, i - left.size());
    return parent.GetNode();
  }

  MersenneTwister mersenne_twister_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("RootwardSampler") {
  Driver driver;
  auto tree_collection = RootedTreeCollection::OfTreeCollection(
  driver.ParseNewickFile("data/five_taxon_rooted.nwk"));
  SubsplitDAG dag(tree_collection);

  RootwardSampler sampler;
  auto tree = sampler.SampleTopology(sampler.RandomLeaf(dag));
  std::cout << tree->Newick([](const Node* node){
    return std::to_string(node->Id());
  }) << "\n";
}

#endif // DOCTEST_LIBRARY_INCLUDED
