// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree_collection.hpp"
#include <regex>
#include "numerical_utils.hpp"

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
  TaxonDateMap taxon_date_map;
  std::regex date_regex(R"raw(^.+_(\d*\.?\d+(?:[eE][-+]?\d+)?)$)raw");
  std::smatch match_date;
  bool first_pass_through_parsing_loop = true;
  bool have_parsed_a_date = false;
  double max_date = DOUBLE_NEG_INF;
  for (auto &[tag, taxon] : TagTaxonMap()) {
    // TODO factor this up
    size_t id = static_cast<size_t>(UnpackFirstInt(tag));
    if (std::regex_match(taxon, match_date, date_regex)) {
      double date = std::stod(match_date[1].str());
      max_date = std::max(date, max_date);
      SafeInsert(taxon_date_map, id, date);
      if (first_pass_through_parsing_loop) {
        have_parsed_a_date = true;
      } else {
        Failwith("We couldn't parse dates for a while, but we could parse:" + taxon);
      }
    } else {  // We couldn't parse a date.
      if (!first_pass_through_parsing_loop && have_parsed_a_date) {
        Failwith("We did parse at least one date, but couldn't parse:" + taxon);
      }  // else this is the first pass through the loop or we haven't parsed a date.
      SafeInsert(taxon_date_map, id, 0.);
    }
    first_pass_through_parsing_loop = false;
  }
  if (have_parsed_a_date) {
    for (auto& [id, date] : taxon_date_map) {
      date -= max_date;
    }
  }
  return taxon_date_map;
}
