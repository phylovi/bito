// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_HPP_
#define SRC_TREE_HPP_

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "node.hpp"
#include "sugar.hpp"

class Tree {
 public:
  typedef std::vector<Tree> TreeVector;
  typedef std::vector<double> BranchLengthVector;

  Tree() {}

  // This is the primary constructor.
  explicit Tree(Node::NodePtr topology, BranchLengthVector branch_lengths);

  // This constructor takes a map of tags to branch lengths; this map gets
  // turned into a branch length vector. It re-ids the topology. Note: any
  // missing branch lengths are set to zero.
  explicit Tree(Node::NodePtr topology, TagDoubleMap branch_lengths);

  const Node::NodePtr Topology() const { return topology_; }
  const BranchLengthVector BranchLengths() const { return branch_lengths_; }
  uint32_t LeafCount() const { return Topology()->LeafCount(); }
  Node::NodePtrVec Children() const { return Topology()->Children(); }
  size_t Id() const { return Topology()->Id(); }
  std::vector<size_t> ParentIdVector() const { return Topology()->ParentIdVector(); }

  bool operator==(const Tree& other) const;

  std::string Newick() const { return Newick(std::nullopt); }
  std::string Newick(TagStringMapOption node_labels) const;

  double BranchLength(const Node* node) const;

  // Take a bifurcating tree and move the root position so that the left hand
  // branch has zero branch length. Modifies tree in place.
  void SlideRootPosition();

  // Returns a new version of this tree without a trifurcation at the root,
  // making it a bifurcation. Given (s0:b0, s1:b1, s2:b2):b4, we get (s0:b0,
  // (s1:b1, s2:b2):0):0. Note that we zero out the root branch length.
  Tree Detrifurcate() const;

  static Tree UnitBranchLengthTreeOf(Node::NodePtr topology);
  static Tree OfParentIdVector(std::vector<size_t> indices);
  static TreeVector ExampleTrees();

  // We make branch lengths public so we can muck with them in Python.
  BranchLengthVector branch_lengths_;

 private:
  Node::NodePtr topology_;
};

inline bool operator!=(const Tree& lhs, const Tree& rhs) { return !(lhs == rhs); }

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Tree") {
  auto trees = Tree::ExampleTrees();
  auto original_newick = trees[0].Newick();
  CHECK_EQ(trees[0].Detrifurcate().Topology(), trees[3].Topology());
  // Shows that Detrifurcate doesn't change the original tree.
  CHECK_EQ(original_newick, trees[0].Newick());
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_TREE_HPP_
