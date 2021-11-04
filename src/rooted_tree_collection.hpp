// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A rooted tree collection has a notion of sampling date for the tips of the tree, and
// all taxa are assumed to share those sampling dates.

#pragma once

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

  void SetDatesToBeConstant(bool initialize_time_trees_using_branch_lengths);
  void ParseDatesFromTaxonNames(bool initialize_time_trees_using_branch_lengths);
  void ParseDatesFromCSV(const std::string& csv_path,
                         bool initialize_time_trees_using_branch_lengths);

 private:
  TagDateMap tag_date_map_;

  void SetTipDates();
  void ProcessTreeDates(bool initialize_time_trees_using_branch_lengths);
  void ParseDatesFromCSVButDontInitializeTimeTrees(const std::string& csv_path);
};

#ifdef DOCTEST_LIBRARY_INCLUDED
// Test of ParseDatesFromTaxonNames appears in rooted_sbn_instance.hpp.
TEST_CASE("RootedTreeCollection") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

