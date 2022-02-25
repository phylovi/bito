// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "unrooted_tree.hpp"

UnrootedTree::UnrootedTree(const Node::NodePtr& topology,
                           BranchLengthVector branch_lengths)
    : Tree(topology, std::move(branch_lengths)) {
  AssertTopologyTrifurcatingInConstructor(topology);
}

UnrootedTree::UnrootedTree(const Node::NodePtr& topology, TagDoubleMap branch_lengths)
    : Tree(topology, std::move(branch_lengths)) {
  AssertTopologyTrifurcatingInConstructor(topology);
}

UnrootedTree UnrootedTree::UnitBranchLengthTreeOf(const Node::NodePtr& topology) {
  topology->Polish();
  return UnrootedTree(topology, BranchLengthVector(1 + topology->Id(), 1.));
}

UnrootedTree UnrootedTree::OfParentIdVector(const std::vector<size_t>& ids) {
  auto topology = Node::OfParentIdVector(ids);
  return UnrootedTree(topology, BranchLengthVector(1 + topology->Id(), 1.));
}

Tree UnrootedTree::Detrifurcate() const {
  Assert(Children().size() == 3,
         "UnrootedTree::Detrifurcate given a non-trifurcating tree.");
  auto branch_lengths = BranchLengths();
  auto our_id = Id();
  auto root12 = Node::Join(Children()[1], Children()[2], our_id);
  branch_lengths[our_id] = 0.;
  auto rerooted_topology = Node::Join(Children()[0], root12, our_id + 1);
  branch_lengths.push_back(0.);
  return Tree(rerooted_topology, branch_lengths);
}

UnrootedTree UnrootedTree::DeepCopy() const {
  return UnrootedTree(Topology()->DeepCopy(), BranchLengths());
}

bool UnrootedTree::operator==(const UnrootedTree& other) const {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}

void UnrootedTree::AssertTopologyTrifurcatingInConstructor(
    const Node::NodePtr& topology) {
  Assert(topology->Children().size() == 3,
         "Expected a tree with a trifucation at the root in the constructor of "
         "UnrootedTree.");
}
