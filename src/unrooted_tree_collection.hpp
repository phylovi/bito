// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "tree_collection.hpp"
#include "unrooted_tree.hpp"

using PreUnrootedTreeCollection = GenericTreeCollection<UnrootedTree>;

class UnrootedTreeCollection : public PreUnrootedTreeCollection {
 public:
  // Inherit all constructors.
  using PreUnrootedTreeCollection::PreUnrootedTreeCollection;
  UnrootedTreeCollection(const PreUnrootedTreeCollection& pre_collection);

  static UnrootedTreeCollection OfTreeCollection(const TreeCollection& trees);
};
