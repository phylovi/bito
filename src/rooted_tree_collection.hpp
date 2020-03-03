// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_TREE_COLLECTION_HPP_
#define SRC_ROOTED_TREE_COLLECTION_HPP_

#include "generic_tree_collection.hpp"
#include "rooted_tree.hpp"
#include "tree_collection.hpp"

template class GenericTreeCollection<RootedTree>;
using PreRootedTreeCollection = GenericTreeCollection<RootedTree>;

class RootedTreeCollection : public PreRootedTreeCollection {
  using TagDateMap = std::unordered_map<Tag, double>;

 public:
  // Inherit all constructors.
  using PreRootedTreeCollection::PreRootedTreeCollection;

  static RootedTreeCollection OfTreeCollection(const TreeCollection& trees);

  void ParseDatesFromTaxonNames();

  TagDateMap tag_date_map_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedTreeCollection") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ROOTED_TREE_COLLECTION_HPP_
