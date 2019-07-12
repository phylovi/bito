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

  explicit TreeCollection(Tree::TreePtrVector trees) : trees_(trees) {
    if (trees.size() > 0) {
      auto leaf_count = trees[0]->LeafCount();
      for (const auto &tree : trees) {
        assert(tree->LeafCount() == leaf_count);
      }
    }
  }
  TreeCollection(Tree::TreePtrVector trees, TagStringMap tag_taxon_map)
      : trees_(trees), tag_taxon_map_(tag_taxon_map) {
    auto taxon_count = tag_taxon_map.size();
    for (const auto &tree : trees) {
      assert(tree->LeafCount() == taxon_count);
    }
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

  Node::TopologyCounter TopologyCounter() {
    Node::TopologyCounter counter;
    for (const auto &tree : trees_) {
      auto search = counter.find(tree->Root());
      if (search == counter.end()) {
        assert(counter.insert(std::make_pair(tree->Root(), 1)).second);
      } else {
        search->second++;
      }
    }
    return counter;
  }

 private:
  Tree::TreePtrVector trees_;
  TagStringMap tag_taxon_map_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TopologyCounter") {
  TreeCollection collection(Tree::ExampleTrees());
  auto counter = collection.TopologyCounter();
  std::vector<uint32_t> v;
  for (const auto &iter : counter) {
    v.push_back(iter.second);
  }
  std::vector<uint32_t> v_correct = {1, 3};
  CHECK_EQ(v, v_correct);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_COLLECTION_HPP_
