// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree.hpp"

RootedTree::RootedTree(const Tree& tree) : Tree(tree.Topology(), tree.BranchLengths()) {
  Assert(Children().size() == 2,
         "Failed to create a RootedTree out of a Tree that isn't bifurcating at the "
         "root. Are you trying to parse unrooted trees into a RootedSBNInstance?");
}

bool RootedTree::operator==(const RootedTree& other) const {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}
