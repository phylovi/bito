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
 public:
  // Inherit all constructors.
  using PreRootedTreeCollection::PreRootedTreeCollection;

  // TODO add a parse_dates flag which populates the taxon_date_map_ variable.
  static RootedTreeCollection OfTreeCollection(const TreeCollection& trees);

 private:
  // TODO add this:
  // std::unordered_map<size_t, double> taxon_date_map_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedTreeCollection") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ROOTED_TREE_COLLECTION_HPP_
