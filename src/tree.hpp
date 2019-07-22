// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_HPP_
#define SRC_TREE_HPP_

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "node.hpp"
#include "typedefs.hpp"

class Tree {
 public:
  typedef std::shared_ptr<Tree> TreePtr;
  typedef std::vector<TreePtr> TreePtrVector;
  typedef std::vector<double> BranchLengthVector;

  // This constructor takes a map of tags to branch lengths; this map gets
  // turned into a branch length vector. It reindexes the topology. Note: any
  // missing branch lengths are set to zero.
  explicit Tree(Node::NodePtr topology, TagDoubleMap branch_lengths);

  explicit Tree(Node::NodePtr topology, BranchLengthVector branch_lengths);

  const Node::NodePtr Topology() const { return topology_; }
  const BranchLengthVector BranchLengths() const { return branch_lengths_; }
  uint32_t LeafCount() const { return Topology()->LeafCount(); }
  Node::NodePtrVec Children() const { return Topology()->Children(); }
  size_t Index() const { return Topology()->Index(); }

  bool operator==(const Tree& other);

  std::string Newick(
      TagStringMapOption node_labels = std::experimental::nullopt) const;

  double BranchLength(const Node* node) const;

  // Take a bifurcating tree and move the root position so that the left hand
  // branch has zero branch length. Modifies tree in place.
  void SlideRootPosition();

  // Returns a new version of this tree without a trifurcation at the root,
  // making it a bifurcation. Given (s0:b0, s1:b1, s2:b2):b4, we get (s0:b0,
  // (s1:b1, s2:b2):0):0. Note that we zero out the root branch length.
  TreePtr Detrifurcate();
  static TreePtr UnitBranchLengthTreeOf(Node::NodePtr topology);
  static TreePtrVector ExampleTrees();

 private:
  Node::NodePtr topology_;
  BranchLengthVector branch_lengths_;
};

// Compare TreePtrs by their Trees.
inline bool operator==(const Tree::TreePtr& lhs, const Tree::TreePtr& rhs) {
  return *lhs == *rhs;
}

inline bool operator!=(const Tree::TreePtr& lhs, const Tree::TreePtr& rhs) {
  return !(lhs == rhs);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Tree") {
  auto trees = Tree::ExampleTrees();
  auto original_newick = trees[0]->Newick();
  CHECK_EQ(trees[0]->Detrifurcate()->Topology(), trees[3]->Topology());
  // Shows that Detrifurcate doesn't change the original tree.
  CHECK_EQ(original_newick, trees[0]->Newick());
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_TREE_HPP_
