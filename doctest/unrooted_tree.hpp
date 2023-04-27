#pragma once

#include "../src/unrooted_tree.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("UnrootedTree") {
  auto trees = Tree::ExampleTrees();
  auto unrooted_tree = UnrootedTree(trees[0]);
  auto original_newick = unrooted_tree.Newick();
  CHECK_EQ(unrooted_tree.Detrifurcate().Topology(), trees[3].Topology());
  // Shows that Detrifurcate doesn't change the original tree.
  CHECK_EQ(original_newick, unrooted_tree.Newick());

  auto topologies = Node::ExampleTopologies();
  // This should work: topology has trifurcation at the root.
  UnrootedTree::UnitBranchLengthTreeOf(topologies[0]);
  // This shouldn't.
  CHECK_THROWS_AS(UnrootedTree::UnitBranchLengthTreeOf(topologies[3]),
                  std::runtime_error&);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
