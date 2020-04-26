// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_UNROOTED_TREE_HPP_
#define SRC_UNROOTED_TREE_HPP_

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

  // Returns a new version of this tree without a trifurcation at the root,
  // making it a bifurcation. Given (s0:b0, s1:b1, s2:b2):b4, we get (s0:b0,
  // (s1:b1, s2:b2):0):0. Note that we zero out the root branch length.
  Tree Detrifurcate() const;

  static UnrootedTree UnitBranchLengthTreeOf(Node::NodePtr topology);
  static UnrootedTree OfParentIdVector(std::vector<size_t> indices);
  static TreeVector ExampleTrees() = delete;

 private:
  static void AssertTopologyTrifurcatingInConstructor(const Node::NodePtr& topology);
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("UnrootedTree") {
  auto trees = Tree::ExampleTrees();
  auto original_newick = trees[0].Newick();
  CHECK_EQ(UnrootedTree(trees[0]).Detrifurcate().Topology(), trees[3].Topology());
  // Shows that Detrifurcate doesn't change the original tree.
  CHECK_EQ(original_newick, trees[0].Newick());

  auto topologies = Node::ExampleTopologies();
  // This works: topology has trifurcation at the root.
  UnrootedTree::UnitBranchLengthTreeOf(topologies[0]);
  // This shouldn't.
  CHECK_THROWS_AS(UnrootedTree::UnitBranchLengthTreeOf(topologies[3]),
                  std::runtime_error&);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_UNROOTED_TREE_HPP_
