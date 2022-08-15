// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tree.hpp"

template <typename TTree>
class GenericTreeCollection {
 protected:
  using TTreeVector = std::vector<TTree>;

 public:
  GenericTreeCollection() = default;

  explicit GenericTreeCollection(TTreeVector trees) : trees_(std::move(trees)) {
    if (!trees_.empty()) {
      auto leaf_count = trees_[0].LeafCount();
      auto different_leaf_count = [leaf_count](const auto &tree) {
        return tree.LeafCount() != leaf_count;
      };
      if (std::any_of(trees_.cbegin(), trees_.cend(), different_leaf_count)) {
        Failwith("Trees must all have the same number of tips in a tree collection.");
      }
    }
  }

  GenericTreeCollection(TTreeVector trees, TagStringMap tag_taxon_map)
      : trees_(std::move(trees)), tag_taxon_map_(std::move(tag_taxon_map)) {
    auto taxon_count = tag_taxon_map.size();
    auto different_taxon_count = [taxon_count](const auto &tree) {
      return tree.LeafCount() != taxon_count;
    };
    if (std::any_of(trees.cbegin(), trees.cend(), different_taxon_count)) {
      Failwith(
          "Tree leaf count doesn't match the size of tag_taxon_map when building a "
          "tree collection.");
    }
  }

  GenericTreeCollection(TTreeVector trees, const std::vector<std::string> &taxon_labels)
      : GenericTreeCollection(std::move(trees), TagStringMapOf(taxon_labels)) {}

  size_t TreeCount() const { return trees_.size(); }
  const TTreeVector &Trees() const { return trees_; }
  const TTree &GetTree(size_t i) const { return trees_.at(i); }
  const TagStringMap &TagTaxonMap() const { return tag_taxon_map_; }
  size_t TaxonCount() const { return tag_taxon_map_.size(); }

  bool operator==(const GenericTreeCollection<TTree> &other) const {
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

  // Remove trees from begin_idx to just before end_idx.
  void Erase(size_t begin_idx, size_t end_idx) {
    if (begin_idx > end_idx || end_idx > TreeCount()) {
      Failwith("Illegal arguments to Tree_Collection.Erase.");
    }
    // else:
    using difference_type = typename TTreeVector::difference_type;
    trees_.erase(trees_.begin() + static_cast<difference_type>(begin_idx),
                 trees_.begin() + static_cast<difference_type>(end_idx));
  }

  // Drop the first fraction trees from the collection.
  void DropFirst(double fraction) {
    Assert(fraction >= 0. && fraction <= 1., "Illegal argument to DropFirst.");
    auto end_idx = static_cast<size_t>(fraction * static_cast<double>(TreeCount()));
    Erase(0, end_idx);
  }

  // Build a tree collection by duplicating the first tree loaded.
  GenericTreeCollection<TTree> BuildCollectionByDuplicatingFirst(
      size_t number_of_times) {
    TTreeVector tree_vector;

    Assert(TreeCount() > 0, "Need at least one tree if we are to duplicate the first.");

    tree_vector.reserve(number_of_times);
    for (size_t idx = 0; idx < number_of_times; idx++) {
      tree_vector.push_back(GetTree(0).DeepCopy());
    }

    return GenericTreeCollection<TTree>(std::move(tree_vector), TagTaxonMap());
  }

  std::string Newick() const {
    std::string str;
    for (const auto &tree : trees_) {
      if (tag_taxon_map_.empty()) {
        str.append(tree.Newick());
      } else {
        str.append(tree.Newick(tag_taxon_map_));
      }
      str.push_back('\n');
    }
    return str;
  }

  void ToNewickFile(const std::string &out_path) const {
    std::ofstream out_stream(out_path);
    out_stream << Newick();
    out_stream.close();
    if (!out_stream) {
      Failwith("ToNewickFile: could not write file to " + out_path);
    }
  }

  void ToNewickTopologyFile(const std::string &out_path) const {
    std::ofstream out_stream(out_path);
    for (const auto &tree : trees_) {
      out_stream << tree.NewickTopology(tag_taxon_map_) << std::endl;
    }
    out_stream.close();
    if (!out_stream) {
      Failwith("ToNewickTopologyFile: could not write file to " + out_path);
    }
  }

  Node::TopologyCounter TopologyCounter() const {
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

  std::vector<std::string> TaxonNames() const {
    std::vector<std::string> names(tag_taxon_map_.size());
    for (const auto &iter : tag_taxon_map_) {
      size_t id = MaxLeafIDOfTag(iter.first);
      Assert(id < names.size(),
             "Leaf ID is out of range in TaxonNames for tree collection.");
      names[id] = iter.second;
    }
    return names;
  }

  static TagStringMap TagStringMapOf(const std::vector<std::string> &taxon_labels) {
    TagStringMap taxon_map;
    for (size_t index = 0; index < taxon_labels.size(); index++) {
      SafeInsert(taxon_map, PackInts(static_cast<uint32_t>(index), 1),
                 taxon_labels[index]);
    }
    return taxon_map;
  }

  static GenericTreeCollection UnitBranchLengthTreesOf(
      std::vector<Node::NodePtr> topologies, TagStringMap tag_taxon_map) {
    std::vector<TTree> tree_vector;
    for (const auto &topology : topologies) {
      tree_vector.push_back(TTree::UnitBranchLengthTreeOf(topology));
    }
    return GenericTreeCollection(tree_vector, tag_taxon_map);
  }

  auto begin() const { return trees_.begin(); }
  auto end() const { return trees_.end(); }

  TTreeVector trees_;

 protected:
  TagStringMap tag_taxon_map_;
};

// Tests appear in non-generic subclasses.
