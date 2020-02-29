// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree_collection.hpp"
#include <regex>

RootedTreeCollection RootedTreeCollection::OfTreeCollection(
    const TreeCollection& trees) {
  TTreeVector rooted_trees;
  rooted_trees.reserve(trees.TreeCount());
  for (const auto& tree : trees.Trees()) {
    rooted_trees.push_back(RootedTree(tree));
  }
  return RootedTreeCollection(std::move(rooted_trees), trees.TagTaxonMap());
}

RootedTreeCollection::TaxonDateMap RootedTreeCollection::ParseDates() {
  std::unordered_map<size_t, double> taxon_date_map;
  std::regex date_regex("^.+_(\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)$");
  std::smatch match_date;
  for (auto &iter : TagTaxonMap()) {
    if (std::regex_match(iter.second, match_date, date_regex)) {
      taxon_date_map.insert(
          std::make_pair(UnpackFirstInt(iter.first), std::stod(match_date[1].str())));
    }
  }
  if (taxon_date_map.size() != 0 && taxon_date_map.size() != TaxonCount()) {
    Failwith("Cannot read dates from tree file.");
  }
  if (taxon_date_map.size() == 0) {
    for (auto &iter : TagTaxonMap()) {
      taxon_date_map.insert(std::make_pair(UnpackFirstInt(iter.first), 0));
    }
  }

  std::vector<double> dates;
  for (const auto &pair : taxon_date_map) {
    dates.push_back(pair.second);
  }
  std::sort(dates.begin(), dates.end());

  // date in years
  if (dates[0] != 0.0) {
    double max = dates[dates.size() - 1];
    for (auto &pair : taxon_date_map) {
      pair.second = max - pair.second;
    }
  }

  return taxon_date_map;
}
