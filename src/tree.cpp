// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

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
    : branch_lengths_(std::move(branch_lengths)), topology_(topology) {
  Assert(topology->Id() + 1 == branch_lengths_.size(),
         "Root id is too large relative to the branch_lengths size in "
         "Tree::Tree.");
}

bool Tree::operator==(const Tree& other) const {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}

std::string Tree::Newick(const TagStringMapOption& node_labels) const {
  return Topology()->Newick(branch_lengths_, node_labels);
}

std::string Tree::NewickTopology(const TagStringMapOption& node_labels) const {
  return Topology()->Newick(std::nullopt, node_labels);
}

double Tree::BranchLength(const Node* node) const {
  Assert(node->Id() < branch_lengths_.size(),
         "Requested id is out of range in Tree::BranchLength.");
  return branch_lengths_[node->Id()];
}

Tree Tree::UnitBranchLengthTreeOf(Node::NodePtr topology) {
  topology->Polish();
  return Tree(topology, BranchLengthVector(1 + topology->Id(), 1.));
}

Tree Tree::OfParentIdVector(const std::vector<size_t>& ids) {
  auto topology = Node::OfParentIdVector(ids);
  return Tree(topology, BranchLengthVector(topology->Id() + 1, 1.));
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
