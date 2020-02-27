// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_TREE_COLLECTION_HPP_
#define SRC_ROOTED_TREE_COLLECTION_HPP_

#include "rooted_tree.hpp"
#include "tree_collection.hpp"

class RootedTreeCollection {
 public:
  RootedTreeCollection() = default;
  explicit RootedTreeCollection(const TreeCollection& trees);

  size_t TreeCount() const { return trees_.size(); }
  const RootedTree& GetTree(size_t i) const { return trees_.at(i); }

  RootedTree::RootedTreeVector trees_;

 private:
  TagStringMap tag_taxon_map_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedTreeCollection") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ROOTED_TREE_COLLECTION_HPP_
