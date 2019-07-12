// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_COLLECTION_HPP_
#define SRC_TREE_COLLECTION_HPP_

#include <memory>
#include <string>
#include "tree.hpp"

class TreeCollection {
 public:
  typedef std::shared_ptr<TreeCollection> TreeCollectionPtr;

  TreeCollection(Tree::TreePtrVector trees, TagStringMap tag_taxon_map)
      : trees_(trees), tag_taxon_map_(tag_taxon_map) {
    // TODO(erick) assert that the number of leaves in the trees is equal to
    // that in tag_taxon_map
  }

  size_t TreeCount() { return trees_.size(); }
  const Tree::TreePtrVector &Trees() const { return trees_; }
  const TagStringMap &TagTaxonMap() const { return tag_taxon_map_; }
  size_t TaxonCount() const { return tag_taxon_map_.size(); }
  std::string Newick() const {
    std::string str;
    for (const auto &tree : trees_) {
      str.append(tree->Newick(tag_taxon_map_));
      str.push_back('\n');
    }
    return str;
  }

 private:
  Tree::TreePtrVector trees_;
  TagStringMap tag_taxon_map_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
// TODO(erick) add tests for tree_collection
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_COLLECTION_HPP_
