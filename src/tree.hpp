// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "node.hpp"
#include "sugar.hpp"

class Tree {
 public:
  using TreeVector = std::vector<Tree>;
  using BranchLengthVector = std::vector<double>;

  Tree() = default;

  // This is the primary constructor. The branch lengths are indexed according to the
  // numbering of the nodes of the tree (see node.hpp for details on how that works.)
  explicit Tree(const Node::NodePtr& topology, BranchLengthVector branch_lengths);

  // This constructor takes a map of tags to branch lengths; this map gets
  // turned into a branch length vector. It re-ids the topology. Note: any
  // missing branch lengths are set to zero.
  explicit Tree(const Node::NodePtr& topology, TagDoubleMap branch_lengths);

  const Node::NodePtr Topology() const { return topology_; }
  const BranchLengthVector& BranchLengths() const { return branch_lengths_; }
  uint32_t LeafCount() const { return Topology()->LeafCount(); }
  Node::NodePtrVec Children() const { return Topology()->Children(); }
  size_t Id() const { return Topology()->Id(); }
  std::vector<size_t> ParentIdVector() const { return Topology()->ParentIdVector(); }

  bool operator==(const Tree& other) const;

  std::string Newick() const { return Newick(std::nullopt); }
  std::string Newick(TagStringMapOption node_labels) const;
  std::string NewickTopology() const;

  double BranchLength(const Node* node) const;

  // Take a bifurcating tree and move the root position so that the left hand
  // branch has zero branch length. Modifies tree in place.
  void SlideRootPosition();

  static Tree UnitBranchLengthTreeOf(Node::NodePtr topology);
  static Tree OfParentIdVector(const std::vector<size_t>& indices);
  static TreeVector ExampleTrees();

  // We make branch lengths public so we can muck with them in Python.
  BranchLengthVector branch_lengths_;

 protected:
  Node::NodePtr topology_;
};

inline bool operator!=(const Tree& lhs, const Tree& rhs) { return !(lhs == rhs); }

#ifdef DOCTEST_LIBRARY_INCLUDED
// Lots of tests in UnrootedTree and RootedTree.
TEST_CASE("Tree") {}
#endif  // DOCTEST_LIBRARY_INCLUDED
