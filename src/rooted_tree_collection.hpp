// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_TREE_COLLECTION_HPP_
#define SRC_ROOTED_TREE_COLLECTION_HPP_

#include "generic_tree_collection.hpp"
#include "rooted_tree.hpp"
#include "tree_collection.hpp"

template class GenericTreeCollection<RootedTree>;
using PreRootedTreeCollection = GenericTreeCollection<RootedTree>;

class RootedTreeCollection : public PreRootedTreeCollection {
  using TagDateMap = TagDoubleMap;

 public:
  // Inherit all constructors.
  using PreRootedTreeCollection::PreRootedTreeCollection;

  static RootedTreeCollection OfTreeCollection(const TreeCollection& trees);

  const TagDateMap& GetTagDateMap() const { return tag_date_map_; };

  void SetDatesToBeConstant(bool initialize_time_trees);
  void ParseDatesFromTaxonNames(bool initialize_time_trees);
  void ParseDatesFromCSV(const std::string& csv_path, bool initialize_time_trees);

 private:
  TagDateMap tag_date_map_;

  void InitializeTimeTrees();
  void SetNodeBoundsUsingDates();
  void PostprocessTrees(bool initialize_time_trees);
  void ParseDatesFromCSVButDontInitializeTimeTrees(const std::string& csv_path);
};

#ifdef DOCTEST_LIBRARY_INCLUDED
// Test of ParseDatesFromTaxonNames appears in rooted_sbn_instance.hpp.
TEST_CASE("RootedTreeCollection") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ROOTED_TREE_COLLECTION_HPP_
