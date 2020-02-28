// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree_collection.hpp"

RootedTreeCollection::RootedTreeCollection(const TreeCollection& trees)
    : tag_taxon_map_(trees.TagTaxonMap()) {
  trees_.reserve(trees.TreeCount());
  for (const auto& tree : trees.Trees()) {
    trees_.push_back(RootedTree(tree));
  }
}

// void SBNInstance::ParseDates() {
//   std::unordered_map<size_t, double> taxon_date_map;
//   std::regex date_regex("^.+_(\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)$");
//   std::smatch match_date;
//   for (auto &iter : tree_collection_.TagTaxonMap()) {
//     if (std::regex_match(iter.second, match_date, date_regex)) {
//       taxon_date_map.insert(
//           std::make_pair(UnpackFirstInt(iter.first),
//           std::stod(match_date[1].str())));
//     }
//   }
//   if (taxon_date_map.size() != 0 &&
//       taxon_date_map.size() != tree_collection_.TaxonCount()) {
//     Failwith("Cannot read dates from tree file.");
//   }
//   if (taxon_date_map.size() == 0) {
//     for (auto &iter : tree_collection_.TagTaxonMap()) {
//       taxon_date_map.insert(std::make_pair(UnpackFirstInt(iter.first), 0));
//     }
//   }
//
//   std::vector<double> dates;
//   for (const auto &pair : taxon_date_map) {
//     dates.push_back(pair.second);
//   }
//   std::sort(dates.begin(), dates.end());
//
//   // date in years
//   if (dates[0] != 0.0) {
//     double max = dates[dates.size() - 1];
//     for (auto &pair : taxon_date_map) {
//       pair.second = max - pair.second;
//     }
//   }
//
//   Tree::TreeVector trees;
//   for (size_t i = 0; i < tree_collection_.trees_.size(); i++) {
//     auto t = std::unique_ptr<Tree>(
//         new RootedTree(*tree_collection_.trees_[i], taxon_date_map));
//     tree_collection_.trees_[i] = std::move(t);
//   }
//   tree_collection_.SetTaxonDateMap(taxon_date_map);
// }

