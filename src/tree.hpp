// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

// TODO(erick)
// make everything const

#ifndef SRC_TREE_HPP_
#define SRC_TREE_HPP_

#include <unordered_map>
#include "node.hpp"

class Tree {
 public:
  typedef std::shared_ptr<Tree> TreePtr;
  typedef std::unordered_map<uint64_t, double> BranchLengthMap;
  typedef std::unordered_map<TreePtr, unsigned int> TreePtrCounter;
  typedef std::shared_ptr<TreePtrCounter> TreePtrCounterPtr;

  explicit Tree(Node::NodePtr root, BranchLengthMap branch_lengths)
      : root_(root), branch_lengths_(branch_lengths) {}

  const Node::NodePtr Root() const { return root_; }
  uint32_t LeafCount() const { return root_->LeafCount(); }
  // TODO(ematsen)
  std::string Newick() const { return root_->Newick(); }

 private:
  Node::NodePtr root_;
  BranchLengthMap branch_lengths_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_HPP_
