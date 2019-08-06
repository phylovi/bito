// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tree.hpp"
#include <cassert>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "node.hpp"
#include "sugar.hpp"

Tree::Tree(Node::NodePtr topology, TagDoubleMap branch_lengths)
    : topology_(topology) {
  auto tag_index_map = topology->Reindex();
  branch_lengths_ = std::vector<double>(topology->Index() + 1);
  for (const auto& iter : tag_index_map) {
    auto& tag = iter.first;
    auto& index = iter.second;
    auto search = branch_lengths.find(tag);
    if (search != branch_lengths.end()) {
      Assert(index < branch_lengths_.size(),
             "branch_lengths of insufficient size in Tree::Tree.");
      branch_lengths_[index] = search->second;
    } else {
      branch_lengths_[index] = 0.;
    }
  }
}

Tree::Tree(Node::NodePtr topology, BranchLengthVector branch_lengths)
    : branch_lengths_(branch_lengths), topology_(topology) {
  Assert(topology->Index() + 1 == branch_lengths.size(),
         "Root index is too large relative to the branch_lengths size in "
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
  Assert(node->Index() < branch_lengths_.size(),
         "Requested index is out of range in Tree::BranchLength.");
  return branch_lengths_[node->Index()];
}

Tree Tree::Detrifurcate() const {
  Assert(Children().size() == 3,
         "Tree::Detrifurcate given a non-trifurcating tree.");
  auto branch_lengths = BranchLengths();
  auto our_index = Index();
  auto root12 = Node::Join(Children()[1], Children()[2], our_index);
  branch_lengths[our_index] = 0.;
  auto rerooted_topology = Node::Join(Children()[0], root12, our_index + 1);
  branch_lengths.push_back(0.);
  return Tree(rerooted_topology, branch_lengths);
}

Tree Tree::UnitBranchLengthTreeOf(Node::NodePtr topology) {
  topology->Reindex();
  BranchLengthVector branch_lengths(1 + topology->Index());
  topology->PreOrder([&branch_lengths](const Node* node) {
    branch_lengths[node->Index()] = 1.;
  });
  return Tree(topology, branch_lengths);
}

Tree Tree::OfParentIndexVector(std::vector<size_t> indices) {
  auto topology = Node::OfParentIndexVector(indices);
  std::vector<double> branch_lengths(topology->Index() + 1, 1.);
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
  size_t fixed_node_index = Children()[1]->Index();
  size_t root_child_index = Children()[0]->Index();
  branch_lengths_[root_child_index] =
      branch_lengths_[root_child_index] + branch_lengths_[fixed_node_index];
  branch_lengths_[fixed_node_index] = 0.0;
}
