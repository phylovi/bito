// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree_collection.hpp"
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

void RootedTreeCollection::ParseDatesFromTaxonNames() {
  tag_date_map_ = TaxonNameMunging::TagDateMapOfTagTaxonMap(TagTaxonMap());
}

void RootedTreeCollection::InitializeParameters() {
  for (auto& tree : trees_) {
    tree.InitializeParameters(tag_date_map_);
  }
}
