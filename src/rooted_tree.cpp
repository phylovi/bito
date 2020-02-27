// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree.hpp"

RootedTree::RootedTree(const Tree& tree)
    : branch_lengths_(tree.BranchLengths()), topology_(tree.Topology()) {}

bool RootedTree::operator==(const RootedTree& other) const {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}

std::string RootedTree::Newick(TagStringMapOption node_labels) const {
  return Topology()->Newick(branch_lengths_, node_labels);
}

double RootedTree::BranchLength(const Node* node) const {
  Assert(node->Id() < branch_lengths_.size(),
         "Requested id is out of range in RootedTree::BranchLength.");
  return branch_lengths_[node->Id()];
}
