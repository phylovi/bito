// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_UNROOTED_TREE_HPP_
#define SRC_UNROOTED_TREE_HPP_

#include <vector>
#include "tree.hpp"

class UnrootedTree : public Tree {
 public:
  typedef std::vector<UnrootedTree> UnrootedTreeVector;

  // See tree.hpp for description of constructors.
  explicit UnrootedTree(const Node::NodePtr& topology,
                        BranchLengthVector branch_lengths);
  explicit UnrootedTree(const Node::NodePtr& topology, TagDoubleMap branch_lengths);

  static UnrootedTree UnitBranchLengthTreeOf(Node::NodePtr topology);

 private:
  static void AssertTopologyTrifurcatingInConstructor(const Node::NodePtr& topology);
};


#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("UnrootedTree") {
  auto topologies = Node::ExampleTopologies();
  // This works: topology has trifurcation at the root.
  UnrootedTree::UnitBranchLengthTreeOf(topologies[0]);
  // This shouldn't.
  CHECK_THROWS_AS(UnrootedTree::UnitBranchLengthTreeOf(topologies[3]),
                  std::runtime_error&);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_UNROOTED_TREE_HPP_
