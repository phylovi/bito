// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tree_collection.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "sugar.hpp"
#include "tree.hpp"

TreeCollection::TreeCollection() {}

TreeCollection::TreeCollection(Tree::TreeVector trees)
    : trees_(std::move(trees)) {
  if (trees.size() > 0) {
    auto leaf_count = trees[0].LeafCount();
    for (const auto &tree : trees) {
      assert(tree.LeafCount() == leaf_count);
    }
  }
}

TreeCollection::TreeCollection(Tree::TreeVector trees,
                               TagStringMap tag_taxon_map)
    : trees_(std::move(trees)), tag_taxon_map_(std::move(tag_taxon_map)) {
  auto taxon_count = tag_taxon_map.size();
  for (const auto &tree : trees) {
    assert(tree.LeafCount() == taxon_count);
  }
}

TreeCollection::TreeCollection(Tree::TreeVector trees,
                               const std::vector<std::string> &taxon_labels)
    : TreeCollection::TreeCollection(
          trees, TreeCollection::TagStringMapOf(taxon_labels)) {}

bool TreeCollection::operator==(const TreeCollection &other) const {
  if (this->TagTaxonMap() != other.TagTaxonMap()) {
    return false;
  }
  if (TreeCount() != other.TreeCount()) {
    return false;
  }
  for (size_t i = 0; i < TreeCount(); i++) {
    if (this->GetTree(i) != other.GetTree(i)) {
      return false;
    }
  }
  return true;
}

std::string TreeCollection::Newick() const {
  std::string str;
  for (const auto &tree : trees_) {
    if (tag_taxon_map_.size()) {
      str.append(tree.Newick(tag_taxon_map_));
    } else {
      str.append(tree.Newick());
    }
    str.push_back('\n');
  }
  return str;
}

Node::TopologyCounter TreeCollection::TopologyCounter() {
  Node::TopologyCounter counter;
  for (const auto &tree : trees_) {
    auto search = counter.find(tree.Topology());
    if (search == counter.end()) {
      SafeInsert(counter, tree.Topology(), static_cast<uint32_t>(1));
    } else {
      search->second++;
    }
  }
  return counter;
}

std::vector<std::string> TreeCollection::TaxonNames() {
  std::vector<std::string> names(tag_taxon_map_.size());
  for (const auto &iter : tag_taxon_map_) {
    size_t id = Node::MaxLeafIDOfTag(iter.first);
    assert(id < names.size());
    names[id] = iter.second;
  }
  return names;
}

TagStringMap TreeCollection::TagStringMapOf(
    std::vector<std::string> taxon_labels) {
  TagStringMap taxon_map;
  for (size_t index = 0; index < taxon_labels.size(); index++) {
    SafeInsert(taxon_map, PackInts(static_cast<uint32_t>(index), 1),
               taxon_labels[index]);
  }
  return taxon_map;
}
