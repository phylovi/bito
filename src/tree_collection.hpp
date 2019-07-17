// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_COLLECTION_HPP_
#define SRC_TREE_COLLECTION_HPP_

#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "tree.hpp"

class TreeCollection {
 public:
  typedef std::shared_ptr<TreeCollection> TreeCollectionPtr;

  explicit TreeCollection(Tree::TreePtrVector trees)
      : trees_(std::move(trees)) {
    if (trees.size() > 0) {
      auto leaf_count = trees[0]->LeafCount();
      for (const auto &tree : trees) {
        assert(tree->LeafCount() == leaf_count);
      }
    }
  }

  TreeCollection(Tree::TreePtrVector trees, TagStringMap tag_taxon_map)
      : trees_(std::move(trees)), tag_taxon_map_(std::move(tag_taxon_map)) {
    auto taxon_count = tag_taxon_map.size();
    for (const auto &tree : trees) {
      assert(tree->LeafCount() == taxon_count);
    }
  }

  size_t TreeCount() const { return trees_.size(); }
  const Tree::TreePtrVector &Trees() const { return trees_; }
  const Tree::TreePtr &GetTree(size_t i) const { return trees_.at(i); }
  const TagStringMap &TagTaxonMap() const { return tag_taxon_map_; }
  size_t TaxonCount() const { return tag_taxon_map_.size(); }

  bool operator==(const TreeCollection &other) {
    if (this->TagTaxonMap() != other.TagTaxonMap()) {
      return false;
    }
    if (TreeCount() != other.TreeCount()) {
      return false;
    }
    for (size_t i = 0; i < TreeCount(); i++) {
      if (this->Trees()[i] != other.Trees()[i]) {
        return false;
      }
    }
    return true;
  }

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
      auto search = counter.find(tree->Topology());
      if (search == counter.end()) {
        assert(counter.insert({tree->Topology(), 1}).second);
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

// Compare TreeCollectionPtrs by their TreeCollections.
inline bool operator==(const TreeCollection::TreeCollectionPtr &lhs,
                       const TreeCollection::TreeCollectionPtr &rhs) {
  return *lhs == *rhs;
}

inline bool operator!=(const TreeCollection::TreeCollectionPtr &lhs,
                       const TreeCollection::TreeCollectionPtr &rhs) {
  return !(lhs == rhs);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TopologyCounter") {
  TreeCollection collection(Tree::ExampleTrees());
  auto counter = collection.TopologyCounter();
  std::unordered_map<std::string, uint32_t> counted;
  for (const auto &iter : counter) {
    assert(counted.insert({iter.first->Newick(), iter.second}).second);
  }
  std::unordered_map<std::string, uint32_t> counted_correct(
      {{"(0_1,1_1,(2_1,3_1)3_2)3_4;", 2},
       {"(0_1,2_1,(1_1,3_1)3_2)3_4;", 1},
       {"(0_1,(1_1,(2_1,3_1)3_2)3_3)3_4;", 1}});
  CHECK_EQ(counted, counted_correct);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_COLLECTION_HPP_
