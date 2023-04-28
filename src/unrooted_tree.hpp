// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <vector>

#include "tree.hpp"

class UnrootedTree : public Tree {
 public:
  typedef std::vector<UnrootedTree> UnrootedTreeVector;

  // See tree.hpp for description of constructors.
  UnrootedTree(const Node::NodePtr& topology, BranchLengthVector branch_lengths);
  UnrootedTree(const Node::NodePtr& topology, TagDoubleMap branch_lengths);
  explicit UnrootedTree(Tree tree)
      : UnrootedTree(tree.Topology(), std::move(tree.branch_lengths_)){};

  UnrootedTree DeepCopy() const;

  bool operator==(const Tree& other) const = delete;
  bool operator==(const UnrootedTree& other) const;

  // Returns a new version of this tree without a trifurcation at the root,
  // making it a bifurcation. Given (s0:b0, s1:b1, s2:b2):b4, we get (s0:b0,
  // (s1:b1, s2:b2):0):0. Note that we zero out the root branch length.
  Tree Detrifurcate() const;

  static UnrootedTree UnitBranchLengthTreeOf(const Node::NodePtr& topology);
  static UnrootedTree OfParentIdVector(const std::vector<size_t>& indices);
  static TreeVector ExampleTrees() = delete;

 private:
  static void AssertTopologyTrifurcatingInConstructor(const Node::NodePtr& topology);
};
