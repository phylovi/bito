// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree.hpp"

RootedTree::RootedTree(const Tree& tree) : Tree(tree.Topology(), tree.BranchLengths()) {
  Assert(Children().size() == 2,
         "Failed to create a RootedTree out of a Tree that isn't bifurcating at the "
         "root. Are you trying to parse unrooted trees into a RootedSBNInstance?");
}

void RootedTree::InitializeParameters(
    const std::unordered_map<Tag, double>& tag_date_map) {
  size_t leaf_count = LeafCount();
  parameters_ = std::vector<double>(leaf_count - 1, -1);

  node_heights_ = std::vector<double>(Topology()->Id() + 1);
  node_bounds_ = std::vector<double>(Topology()->Id() + 1);

  std::vector<double> branch_lengths = BranchLengths();
  int root_id = static_cast<int>(Topology()->Id());

  // Set up bounds and heights
  // Start with the leaves from the dates map
  for (const auto& [tag, date] : tag_date_map) {
    size_t id = MaxLeafIDOfTag(tag);
    node_heights_[id] = date;
    node_bounds_[id] = date;
  }

  // Set up the internal heights and bounds
  Topology()->BinaryIdPostOrder(
      [&leaf_count, &branch_lengths, this](int node_id, int child0_id, int child1_id) {
        if (node_id >= leaf_count) {
          node_bounds_[node_id] =
              std::max(node_bounds_[child0_id], node_bounds_[child1_id]);
          node_heights_[node_id] = node_heights_[child0_id] + branch_lengths[child0_id];
        }
      });

  // Set up ratios
  parameters_[root_id - leaf_count] = node_heights_[root_id];
  Topology()->TripleIdPreOrderBifurcating(
      [&root_id, &leaf_count, this](int node_id, int sister_id, int parent_id) {
        if (node_id >= leaf_count) {
          parameters_[node_id - leaf_count] =
              (node_heights_[node_id] - node_bounds_[node_id]) /
              (node_heights_[parent_id] - node_bounds_[node_id]);
        }
      });
}

bool RootedTree::operator==(const RootedTree& other) const {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}
