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

