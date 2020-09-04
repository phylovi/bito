// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree_collection.hpp"
#include "csv.h"
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
  PostprocessTrees(initialize_time_trees);
}

void RootedTreeCollection::ParseDatesFromTaxonNames(bool initialize_time_trees) {
  tag_date_map_ = TaxonNameMunging::ParseDatesFromTagTaxonMap(TagTaxonMap());
  PostprocessTrees(initialize_time_trees);
}

void RootedTreeCollection::ParseDatesFromCSV(const std::string& csv_path,
                                             bool initialize_time_trees) {
  ParseDatesFromCSVButDontInitializeTimeTrees(csv_path);
  PostprocessTrees(initialize_time_trees);
}

// Given a headerless 2-column CSV of quoted string keys then double values, this parses
// the CSV into a StringDoubleMap.
StringDoubleMap StringDoubleMapOfCSV(const std::string& csv_path) {
  io::CSVReader<2, io::trim_chars<' ', '\t'>, io::double_quote_escape<',', '"'>> csv_in(
      csv_path);
  std::string key;
  double value;
  StringDoubleMap string_double_map;
  while (csv_in.read_row(key, value)) {
    auto search = string_double_map.find(key);
    if (search == string_double_map.end()) {
      string_double_map.insert({key, value});
    } else {
      Failwith("Key " + key + " found twice in " + csv_path +  // NOLINT
               "when turning it into a map.");
    }
  }
  return string_double_map;
}

void RootedTreeCollection::ParseDatesFromCSVButDontInitializeTimeTrees(
    const std::string& csv_path) {
  tag_date_map_.clear();
  auto taxon_date_map = StringDoubleMapOfCSV(csv_path);
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

void RootedTreeCollection::PostprocessTrees(bool initialize_time_trees) {
  if (initialize_time_trees) {
    InitializeTimeTrees();
  } else {
    SetNodeBoundsUsingDates();
  }
}
