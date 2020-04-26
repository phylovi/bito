// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "unrooted_tree.hpp"

UnrootedTree::UnrootedTree(const Node::NodePtr& topology,
                           BranchLengthVector branch_lengths)
    : Tree(topology, branch_lengths) {
  AssertTopologyTrifurcatingInConstructor(topology);
}

UnrootedTree::UnrootedTree(const Node::NodePtr& topology, TagDoubleMap branch_lengths)
    : Tree(topology, branch_lengths) {
  AssertTopologyTrifurcatingInConstructor(topology);
}

// TODO what is up with this tree traversal?
UnrootedTree UnrootedTree::UnitBranchLengthTreeOf(Node::NodePtr topology) {
  topology->Polish();
  BranchLengthVector branch_lengths(1 + topology->Id());
  topology->PreOrder(
      [&branch_lengths](const Node* node) { branch_lengths[node->Id()] = 1.; });
  return UnrootedTree(topology, branch_lengths);
}

UnrootedTree UnrootedTree::OfParentIdVector(std::vector<size_t> ids) {
  auto topology = Node::OfParentIdVector(ids);
  std::vector<double> branch_lengths(topology->Id() + 1, 1.);
  return UnrootedTree(topology, std::move(branch_lengths));
}

void UnrootedTree::AssertTopologyTrifurcatingInConstructor(
    const Node::NodePtr& topology) {
  Assert(topology->Children().size() == 3,
         "Expected a tree with a trifucation at the root in the constructor of "
         "UnrootedTree.");
}
