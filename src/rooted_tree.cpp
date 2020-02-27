// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree.hpp"

RootedTree::RootedTree(const Tree& tree)
    : Tree(tree.Topology(), tree.BranchLengths()) {}

bool RootedTree::operator==(const RootedTree& other) const {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}
