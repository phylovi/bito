// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "generic_tree_collection.hpp"
#include "tree.hpp"

using TreeCollection = GenericTreeCollection<Tree>;

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TreeCollection") {
  auto first_four_example_trees = Tree::ExampleTrees();
  first_four_example_trees.resize(4);
  TreeCollection collection(first_four_example_trees);
  auto counter = collection.TopologyCounter();
  std::unordered_map<std::string, uint32_t> counted;
  for (const auto &iter : counter) {
    SafeInsert(counted, iter.first->Newick(std::nullopt, std::nullopt, true),
               iter.second);
  }
  std::unordered_map<std::string, uint32_t> counted_correct(
      {{"(0_1,1_1,(2_1,3_1)3_2)3_4;", 2},
       {"(0_1,2_1,(1_1,3_1)3_2)3_4;", 1},
       {"(0_1,(1_1,(2_1,3_1)3_2)3_3)3_4;", 1}});
  CHECK_EQ(counted, counted_correct);
  collection.DropFirst(0.25);
  CHECK_EQ(collection.TreeCount(), 3);
  collection.DropFirst(1.);
  CHECK_EQ(collection.TreeCount(), 0);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
