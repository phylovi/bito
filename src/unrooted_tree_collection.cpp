// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "unrooted_tree_collection.hpp"

// Explicit class template instantiation:
// https://en.cppreference.com/w/cpp/language/class_template#Explicit_instantiation
template class GenericTreeCollection<UnrootedTree>;

UnrootedTreeCollection UnrootedTreeCollection::OfTreeCollection(
    const TreeCollection& trees) {
  TTreeVector unrooted_trees;
  unrooted_trees.reserve(trees.TreeCount());
  for (const auto& tree : trees.Trees()) {
    unrooted_trees.emplace_back(tree);
  }
  return UnrootedTreeCollection(std::move(unrooted_trees), trees.TagTaxonMap());
}
