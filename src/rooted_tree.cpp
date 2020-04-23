// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree.hpp"

// Branch lengths should correspond to a time tree up to this tolerance.
constexpr double BRANCH_LENGTH_TOLERANCE = 1e-6;

RootedTree::RootedTree(const Tree& tree) : Tree(tree.Topology(), tree.BranchLengths()) {
  Assert(Children().size() == 2,
         "Failed to create a RootedTree out of a Tree that isn't bifurcating at the "
         "root. Are you trying to parse unrooted trees into a RootedSBNInstance?");
}

void RootedTree::InitializeParameters(
    const std::unordered_map<Tag, double>& tag_date_map) {
  size_t leaf_count = LeafCount();
  int root_id = static_cast<int>(Topology()->Id());
  height_ratios_ = std::vector<double>(leaf_count - 1, -1);
  node_heights_ = std::vector<double>(Topology()->Id() + 1);
  node_bounds_ = std::vector<double>(Topology()->Id() + 1);

  // First initialize the leaves using the date map.
  for (const auto& [tag, date] : tag_date_map) {
    size_t leaf_id = MaxLeafIDOfTag(tag);
    node_heights_[leaf_id] = date;
    node_bounds_[leaf_id] = date;
  }

  // Initialize the internal heights and bounds.
  Topology()->BinaryIdPostOrder(
      [&leaf_count, this](int node_id, int child0_id, int child1_id) {
        if (node_id >= leaf_count) {
          node_bounds_[node_id] =
              std::max(node_bounds_[child0_id], node_bounds_[child1_id]);
          node_heights_[node_id] =
              node_heights_[child0_id] + branch_lengths_[child0_id];
          if (fabs(node_heights_[child1_id] + branch_lengths_[child1_id] -
                   node_heights_[node_id]) > BRANCH_LENGTH_TOLERANCE) {
            Failwith("Tree not ultrametric in RootedTree::InitializeParameters.");
          }
        }
      });

  // Initialize ratios.
  // The "height ratio" for the root is the root height.
  height_ratios_[root_id - leaf_count] = node_heights_[root_id];
  Topology()->TripleIdPreOrderBifurcating(
      [&leaf_count, this](int node_id, int, int parent_id) {
        if (node_id >= leaf_count) {
          height_ratios_[node_id - leaf_count] =
              (node_heights_[node_id] - node_bounds_[node_id]) /
              (node_heights_[parent_id] - node_bounds_[node_id]);
        }
      });
}

TagDoubleMap RootedTree::TagDateMapOfDateVector(std::vector<double> leaf_date_vector) {
  Assert(leaf_date_vector.size() == LeafCount(),
         "Wrong size vector in TagDateMapOfDateVector");
  TagDoubleMap tag_date_map;
  for (uint32_t leaf_id = 0; leaf_id < LeafCount(); ++leaf_id) {
    SafeInsert(tag_date_map, PackInts(leaf_id, 1), leaf_date_vector[leaf_id]);
  }
  return tag_date_map;
}

bool RootedTree::operator==(const RootedTree& other) const {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}
