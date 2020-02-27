// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_TREE_HPP_
#define SRC_ROOTED_TREE_HPP_

#include "tree.hpp"

class RootedTree {
 public:
  typedef std::vector<RootedTree> RootedTreeVector;

  explicit RootedTree(const Tree& tree);

  const Node::NodePtr Topology() const { return topology_; }
  const Tree::BranchLengthVector BranchLengths() const { return branch_lengths_; }
  uint32_t LeafCount() const { return Topology()->LeafCount(); }
  Node::NodePtrVec Children() const { return Topology()->Children(); }
  size_t Id() const { return Topology()->Id(); }
  std::vector<size_t> ParentIdVector() const { return Topology()->ParentIdVector(); }

  bool operator==(const RootedTree& other) const;

  std::string Newick() const { return Newick(std::nullopt); }
  std::string Newick(TagStringMapOption node_labels) const;

  double BranchLength(const Node* node) const;

  // We make branch lengths public so we can muck with them in Python.
  Tree::BranchLengthVector branch_lengths_;

 private:
  Node::NodePtr topology_;
};

inline bool operator!=(const RootedTree& lhs, const RootedTree& rhs) {
  return !(lhs == rhs);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedTree") {
  // TODO
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_ROOTED_TREE_HPP_
