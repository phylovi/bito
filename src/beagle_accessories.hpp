// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// BeagleAccessories are collections of artifacts that we can make in constant
// time given the tree, and remain const througout any operation-gathering tree
// traversal.

#ifndef SRC_BEAGLE_ACCESSORIES_HPP_
#define SRC_BEAGLE_ACCESSORIES_HPP_

#include <numeric>
#include <vector>
#include "libhmsbeagle/beagle.h"
#include "tree.hpp"

struct BeagleAccessories {
  const int beagle_instance_;
  const bool rescaling_;
  const int fixed_node_id_;
  const int root_child_id_;
  const int node_count_;
  const int taxon_count_;
  const int internal_count_;
  // We're not exacty sure what this mysterious_count argument is for.
  // The beagle docs say: Number of partialsBuffer to integrate (input)
  // In the BEASTs it's hardcoded to 1 and in MrBayes it appears to be for
  // covarion models.
  const int mysterious_count_ = 1;
  // Scaling factors are recomputed every time so we don't read them
  // using destinationScaleRead.
  const int destinationScaleRead_ = BEAGLE_OP_NONE;
  // This is the entry of scaleBuffer in which we store accumulated factors.
  const std::vector<int> cumulative_scale_index_;
  const std::vector<int> node_indices_;
  // pattern weights
  const std::vector<int> category_weight_index_ = {0};
  // state frequencies
  const std::vector<int> state_frequency_index_ = {0};
  // indices of parent partialsBuffers
  std::vector<int> upper_partials_index_ = {0};
  // indices of child partialsBuffers
  std::vector<int> node_partial_indices_ = {0};
  // transition probability matrices
  std::vector<int> node_mat_indices_ = {0};
  // first derivative matrices
  std::vector<int> node_deriv_index_ = {0};

  BeagleAccessories(int beagle_instance, bool rescaling, const Tree &tree)
      : beagle_instance_(beagle_instance),
        rescaling_(rescaling),
        fixed_node_id_(static_cast<int>(tree.Topology()->Children()[1]->Id())),
        root_child_id_(static_cast<int>(tree.Topology()->Children()[0]->Id())),
        node_count_(static_cast<int>(tree.BranchLengths().size())),
        taxon_count_(static_cast<int>(tree.LeafCount())),
        internal_count_(taxon_count_ - 1),
        cumulative_scale_index_({rescaling ? 0 : BEAGLE_OP_NONE}),
        node_indices_(IotaVector(node_count_ - 1, 0)) {}

  static std::vector<int> IotaVector(size_t size, int start_value) {
    std::vector<int> v(size);
    std::iota(v.begin(), v.end(), start_value);
    return v;
  }
};

#endif  // SRC_BEAGLE_ACCESSORIES_HPP_
