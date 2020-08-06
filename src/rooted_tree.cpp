// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_tree.hpp"

// Branch lengths should correspond to a time tree up to this tolerance.
constexpr double BRANCH_LENGTH_TOLERANCE = 1e-4;

RootedTree::RootedTree(const Node::NodePtr& topology, BranchLengthVector branch_lengths)
    : Tree(topology, std::move(branch_lengths)) {
  Assert(
      Children().size() == 2,
      "Failed to create a RootedTree out of a topology that isn't bifurcating at the "
      "root. Perhaps you are trying to parse unrooted trees into a RootedSBNInstance?");
}

RootedTree::RootedTree(const Tree& tree)
    : RootedTree(tree.Topology(), tree.BranchLengths()) {}

void RootedTree::InitializeParameters(
    const std::unordered_map<Tag, double>& tag_date_map) {
  size_t leaf_count = LeafCount();
  int root_id = static_cast<int>(Topology()->Id());
  height_ratios_ = std::vector<double>(leaf_count - 1, -1);
  node_heights_ = std::vector<double>(Topology()->Id() + 1);
  node_bounds_ = std::vector<double>(Topology()->Id() + 1);
  rates_ = std::vector<double>(Topology()->Id(), 1.0);
  rate_count_ = 1.0;  // Default is a strict clock with rate 1

  // First initialize the leaves using the date map.
  for (const auto& [tag, date] : tag_date_map) {
    size_t leaf_id = MaxLeafIDOfTag(tag);
    node_heights_[leaf_id] = date;
    node_bounds_[leaf_id] = date;
  }

  // Initialize the internal heights and bounds.
  Topology()->BinaryIdPostOrder([&leaf_count, this](int node_id, int child0_id,
                                                    int child1_id) {
    if (node_id >= leaf_count) {
      node_bounds_[node_id] =
          std::max(node_bounds_[child0_id], node_bounds_[child1_id]);
      node_heights_[node_id] = node_heights_[child0_id] + branch_lengths_[child0_id];
      const auto height_difference =
          fabs(node_heights_[child1_id] + branch_lengths_[child1_id] -
               node_heights_[node_id]);
      if (height_difference > BRANCH_LENGTH_TOLERANCE) {
        Failwith(
            "Tree isn't time-calibrated in RootedTree::InitializeParameters. Height "
            "difference: " +
            std::to_string(height_difference));
      }
    }
  });

  // Initialize ratios.
  // The "height ratio" for the root is the root height.
  height_ratios_[root_id - leaf_count] = node_heights_[root_id];
  Topology()->TripleIdPreOrderBifurcating(
      [&leaf_count, this](int node_id, int, int parent_id) {
        if (node_id >= leaf_count) {
          // See the beginning of the header file for an explanation.
          height_ratios_[node_id - leaf_count] =
              (node_heights_[node_id] - node_bounds_[node_id]) /
              (node_heights_[parent_id] - node_bounds_[node_id]);
        }
      });
}

void RootedTree::SetNodeHeightsViaHeightRatios(EigenConstVectorXdRef height_ratios) {
  size_t leaf_count = LeafCount();
  size_t root_id = Topology()->Id();
  node_heights_[root_id] = height_ratios(root_id - leaf_count);
  Topology()->TripleIdPreOrderBifurcating([&leaf_count, &height_ratios, this](
                                              int node_id, int, int parent_id) {
    if (node_id >= leaf_count) {
      node_heights_[node_id] = node_bounds_[node_id] +
                               height_ratios(node_id - leaf_count) *
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

RootedTree RootedTree::Example() {
  auto topology = Node::ExampleTopologies()[3];
  RootedTree tree(Tree(topology, {2., 1.5, 2., 1., 2.5, 2.5, 0.}));
  std::vector<double> date_vector({5., 3., 0., 1.});
  auto tag_date_map = tree.TagDateMapOfDateVector(date_vector);
  tree.InitializeParameters(tag_date_map);
  return tree;
}

bool RootedTree::operator==(const RootedTree& other) const {
  return (this->Topology() == other.Topology()) &&
         (this->BranchLengths() == other.BranchLengths());
}
