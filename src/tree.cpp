// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tree.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "node.hpp"
#include "sugar.hpp"

Tree::Tree(const Node::NodePtr& topology, TagDoubleMap branch_lengths)
    : topology_(topology) {
  auto tag_id_map = topology->Polish();
  branch_lengths_ = std::vector<double>(topology->Id() + 1);
  for (const auto& [tag, id] : tag_id_map) {
    auto search = branch_lengths.find(tag);
    if (search != branch_lengths.end()) {
      Assert(id < branch_lengths_.size(),
             "branch_lengths of insufficient size in Tree::Tree.");
      branch_lengths_[id] = search->second;
    } else {
      branch_lengths_[id] = 0.;
    }
  }
}

Tree::Tree(const Node::NodePtr& topology, BranchLengthVector branch_lengths)
    : branch_lengths_(branch_lengths), topology_(topology) {
  Assert(topology->Id() + 1 == branch_lengths.size(),
         "Root id is too large relative to the branch_lengths size in "
         "Tree::Tree.");
}

bool Tree::operator==(const Tree& other) const {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}

std::string Tree::Newick(TagStringMapOption node_labels) const {
  return Topology()->Newick(branch_lengths_, node_labels);
}

double Tree::BranchLength(const Node* node) const {
  Assert(node->Id() < branch_lengths_.size(),
         "Requested id is out of range in Tree::BranchLength.");
  return branch_lengths_[node->Id()];
}

Tree Tree::Detrifurcate() const {
  Assert(Children().size() == 3, "Tree::Detrifurcate given a non-trifurcating tree.");
  auto branch_lengths = BranchLengths();
  auto our_id = Id();
  auto root12 = Node::Join(Children()[1], Children()[2], our_id);
  branch_lengths[our_id] = 0.;
  auto rerooted_topology = Node::Join(Children()[0], root12, our_id + 1);
  branch_lengths.push_back(0.);
  return Tree(rerooted_topology, branch_lengths);
}

Tree Tree::UnitBranchLengthTreeOf(Node::NodePtr topology) {
  topology->Polish();
  BranchLengthVector branch_lengths(1 + topology->Id());
  topology->PreOrder(
      [&branch_lengths](const Node* node) { branch_lengths[node->Id()] = 1.; });
  return Tree(topology, branch_lengths);
}

Tree Tree::OfParentIdVector(std::vector<size_t> ids) {
  auto topology = Node::OfParentIdVector(ids);
  std::vector<double> branch_lengths(topology->Id() + 1, 1.);
  return Tree(topology, std::move(branch_lengths));
}

Tree::TreeVector Tree::ExampleTrees() {
  TreeVector v;
  for (const auto& topology : Node::ExampleTopologies()) {
    v.push_back(UnitBranchLengthTreeOf(topology));
  }
  return v;
}

void Tree::SlideRootPosition() {
  size_t fixed_node_id = Children()[1]->Id();
  size_t root_child_id = Children()[0]->Id();
  branch_lengths_[root_child_id] =
      branch_lengths_[root_child_id] + branch_lengths_[fixed_node_id];
  branch_lengths_[fixed_node_id] = 0.0;
}
