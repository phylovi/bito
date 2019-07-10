// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_COLLECTION_HPP_
#define SRC_TREE_COLLECTION_HPP_

#include <memory>
#include <string>
#include <unordered_map>
#include "tree.hpp"

class TreeCollection {
 public:
  typedef std::shared_ptr<TreeCollection> TreeCollectionPtr;
  typedef std::unordered_map<Tree::TreePtr, uint32_t> TreePtrCounter;
  typedef std::shared_ptr<TreePtrCounter> TreePtrCounterPtr;

  explicit TreeCollection(TreePtrCounterPtr trees, TagStringMap tag_taxon_map)
      : trees_(trees), tag_taxon_map_(tag_taxon_map) {}

  size_t TreeCount() { return trees_->size(); }
  TreePtrCounterPtr Trees() const { return trees_; }
  TagStringMap TagTaxonMap() const { return tag_taxon_map_; }
  Tree::TreePtr FirstTree() const {
    assert(trees_->size() > 0);
    return trees_->begin()->first;
  }

 private:
  TreePtrCounterPtr trees_;
  TagStringMap tag_taxon_map_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
// TODO(erick) add tests for tree_collection
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_COLLECTION_HPP_
