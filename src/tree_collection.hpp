// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_COLLECTION_HPP_
#define SRC_TREE_COLLECTION_HPP_

#include "generic_tree_collection.hpp"
#include "tree.hpp"

using TreeCollection = GenericTreeCollection<Tree>;

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TreeCollection") {
  TreeCollection collection(Tree::ExampleTrees());
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

#endif  // SRC_TREE_COLLECTION_HPP_
