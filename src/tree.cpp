// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tree.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "node.hpp"
#include "typedefs.hpp"

Tree::Tree(Node::NodePtr topology, TagDoubleMap branch_lengths)
    : topology_(topology) {
  auto tag_index_map = topology->Reindex();
  branch_lengths_ = std::vector<double>(topology->Index() + 1);
  for (const auto& iter : tag_index_map) {
    auto& tag = iter.first;
    auto& index = iter.second;
    auto search = branch_lengths.find(tag);
    if (search != branch_lengths.end()) {
      assert(index < branch_lengths_.size());
      branch_lengths_[index] = search->second;
    } else {
      branch_lengths_[index] = 0.;
    }
  }
}

Tree::Tree(Node::NodePtr topology, BranchLengthVector branch_lengths)
    : topology_(topology), branch_lengths_(branch_lengths) {
  assert(topology->Index() + 1 == branch_lengths.size());
}

bool Tree::operator==(const Tree& other) {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}

std::string Tree::Newick(TagStringMapOption node_labels) const {
  return Topology()->Newick(branch_lengths_, node_labels);
}

double Tree::BranchLength(const Node* node) const {
  assert(node->Index() < branch_lengths_.size());
  return branch_lengths_[node->Index()];
}

Tree::TreePtr Tree::Detrifurcate() {
  if (Children().size() != 3) {
    std::cerr << "Detrifurcate given a non-trifurcating tree.\n";
    abort();
  }  // else
  auto branch_lengths = BranchLengths();
  auto our_index = Index();
  auto root12 = Node::Join(Children()[1], Children()[2], our_index);
  branch_lengths[our_index] = 0.;
  auto rerooted_topology = Node::Join(Children()[0], root12, our_index + 1);
  branch_lengths.push_back(0.);
  return std::make_shared<Tree>(rerooted_topology, branch_lengths);
}

Tree::TreePtr Tree::UnitBranchLengthTreeOf(Node::NodePtr topology) {
  topology->Reindex();
  BranchLengthVector branch_lengths(1 + topology->Index());
  topology->PreOrder([&branch_lengths](const Node* node) {
    branch_lengths[node->Index()] = 1.;
  });
  return std::make_shared<Tree>(topology, branch_lengths);
}

Tree::TreePtrVector Tree::ExampleTrees() {
  TreePtrVector v;
  for (const auto& topology : Node::ExampleTopologies()) {
    v.push_back(UnitBranchLengthTreeOf(topology));
  }
  return v;
}

// For unrooted, the node with index 1 below the root should have 0 branch
// length. if not
// its value is added to its sibling and set to 0
void Tree::SlideRootPosition() {
  size_t fixed_node_index = Children()[1]->Index();
  size_t root_child_index = Children()[0]->Index();
  branch_lengths_[root_child_index] =
      branch_lengths_[root_child_index] + branch_lengths_[fixed_node_index];
  branch_lengths_[fixed_node_index] = 0.0;
}
