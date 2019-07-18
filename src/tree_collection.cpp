// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tree_collection.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "tree.hpp"

TreeCollection::TreeCollection(Tree::TreePtrVector trees)
    : trees_(std::move(trees)) {
  if (trees.size() > 0) {
    auto leaf_count = trees[0]->LeafCount();
    for (const auto &tree : trees) {
      assert(tree->LeafCount() == leaf_count);
    }
  }
}

TreeCollection::TreeCollection(Tree::TreePtrVector trees,
                               TagStringMap tag_taxon_map)
    : trees_(std::move(trees)), tag_taxon_map_(std::move(tag_taxon_map)) {
  auto taxon_count = tag_taxon_map.size();
  for (const auto &tree : trees) {
    assert(tree->LeafCount() == taxon_count);
  }
}

bool TreeCollection::operator==(const TreeCollection &other) {
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

std::string TreeCollection::Newick() const {
  std::string str;
  for (const auto &tree : trees_) {
    str.append(tree->Newick(tag_taxon_map_));
    str.push_back('\n');
  }
  return str;
}

Node::TopologyCounter TreeCollection::TopologyCounter() {
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
