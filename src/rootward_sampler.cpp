// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "rootward_sampler.hpp"

Node::NodePtr RootwardSampler::SampleTopology(
    SubsplitDAGNode node, EigenConstVectorXdRef inverted_probabilities) {
  Assert(node.IsLeaf(), "Non-leaf node for rootward sampling");
  std::vector<size_t> edges;
  SampleTopology(node, edges, inverted_probabilities);
  Node::NodePtr tree = Node::Leaf(node.Id());
  for (auto i : edges) {
    tree = Node::Join({tree}, i);
  }
  return tree;
}

void RootwardSampler::SetSeed(uint64_t seed) { mersenne_twister_.SetSeed(seed); }

void RootwardSampler::SampleTopology(SubsplitDAGNode node, std::vector<size_t>& edges,
                                     EigenConstVectorXdRef inverted_probabilities) {
  auto left = node.GetLeftRootward();
  auto right = node.GetRightRootward();
  if (left.empty() && right.empty()) {
    // reached root
    return;
  }
  auto parent = SampleParent(left, right, inverted_probabilities);
  edges.push_back(parent.Id());
  SampleTopology(parent, edges, inverted_probabilities);
}

SubsplitDAGNode RootwardSampler::SampleParent(
    ConstNeighborsView left, ConstNeighborsView right,
    EigenConstVectorXdRef inverted_probabilities) {
  std::vector<double> weights;
  weights.resize(left.size() + right.size());
  size_t i = 0;
  for (auto parent = left.begin(); parent != left.end(); ++parent)
    weights[i++] = inverted_probabilities[parent.GetEdgeId()];
  for (auto parent = right.begin(); parent != right.end(); ++parent)
    weights[i++] = inverted_probabilities[parent.GetEdgeId()];
  std::discrete_distribution<> distribution(weights.begin(), weights.end());
  i = static_cast<size_t>(distribution(mersenne_twister_.GetGenerator()));
  if (i < left.size()) {
    auto parent = left.begin();
    std::advance(parent, i);
    return parent.GetNode();
  }
  auto parent = right.begin();
  std::advance(parent, i - left.size());
  return parent.GetNode();
}
