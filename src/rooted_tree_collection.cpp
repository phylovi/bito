// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree_collection.hpp"
#include "csv.hpp"
#include "taxon_name_munging.hpp"

RootedTreeCollection RootedTreeCollection::OfTreeCollection(
    const TreeCollection& trees) {
  TTreeVector rooted_trees;
  rooted_trees.reserve(trees.TreeCount());
  for (const auto& tree : trees.Trees()) {
    rooted_trees.emplace_back(tree);
  }
  return RootedTreeCollection(std::move(rooted_trees), trees.TagTaxonMap());
}

void RootedTreeCollection::SetDatesToBeConstant(bool initialize_time_trees) {
  tag_date_map_ = TaxonNameMunging::ConstantDatesForTagTaxonMap(TagTaxonMap());
  ProcessTreeDates(initialize_time_trees);
}

void RootedTreeCollection::ParseDatesFromTaxonNames(bool initialize_time_trees) {
  tag_date_map_ = TaxonNameMunging::ParseDatesFromTagTaxonMap(TagTaxonMap());
  ProcessTreeDates(initialize_time_trees);
}

void RootedTreeCollection::ParseDatesFromCSV(const std::string& csv_path,
                                             bool initialize_time_trees) {
  ParseDatesFromCSVButDontInitializeTimeTrees(csv_path);
  ProcessTreeDates(initialize_time_trees);
}

void RootedTreeCollection::ParseDatesFromCSVButDontInitializeTimeTrees(
    const std::string& csv_path) {
  tag_date_map_.clear();
  auto taxon_date_map = CSV::StringDoubleMapOfCSV(csv_path);
  for (auto& [tag, taxon] : TagTaxonMap()) {
    auto search = taxon_date_map.find(taxon);
    if (search == taxon_date_map.end()) {
      Failwith("Taxon " + taxon +  // NOLINT
               " found in current tree collection but not in " + csv_path);
    }
    // else
    SafeInsert(tag_date_map_, tag, search->second);
  }
  TaxonNameMunging::MakeDatesRelativeToMaximum(tag_date_map_);
}

void RootedTreeCollection::InitializeTimeTrees() {
  for (auto& tree : trees_) {
    tree.InitializeTimeTree(tag_date_map_);
  }
}

void RootedTreeCollection::SetNodeBoundsUsingDates() {
  for (auto& tree : trees_) {
    tree.SetNodeBoundsUsingDates(tag_date_map_);
  }
}

void RootedTreeCollection::ProcessTreeDates(bool initialize_time_trees) {
  if (initialize_time_trees) {
    InitializeTimeTrees();
  } else {
    SetNodeBoundsUsingDates();
  }
}
