// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tidy_subsplit_dag.hpp"

TidySubsplitDAG::TidySubsplitDAG() : SubsplitDAG() {}

TidySubsplitDAG::TidySubsplitDAG(const RootedTreeCollection &tree_collection)
    : SubsplitDAG(tree_collection) {}

TidySubsplitDAG::TidySubsplitDAG(size_t node_count)
    : above_sorted_(EigenMatrixXb::Identity(node_count, node_count)),
      above_rotated_(EigenMatrixXb::Identity(node_count, node_count)){};

EigenArrayXb TidySubsplitDAG::BelowNode(size_t node_idx) {
  return BelowNode(false, node_idx).max(BelowNode(true, node_idx));
}

EigenArrayXbRef TidySubsplitDAG::BelowNode(bool rotated, size_t node_idx) {
  if (rotated) {
    return above_rotated_.col(node_idx).array();
  } else {
    return above_sorted_.col(node_idx).array();
  }
}

EigenArrayXb TidySubsplitDAG::AboveNode(size_t node_idx) {
  return AboveNode(false, node_idx).max(AboveNode(true, node_idx));
}

EigenArrayXb TidySubsplitDAG::AboveNode(bool rotated, size_t node_idx) {
  if (rotated) {
    return above_rotated_.row(node_idx).array();
  } else {
    return above_sorted_.row(node_idx).array();
  }
}

void TidySubsplitDAG::SetBelow(size_t dst_idx, bool dst_rotated, size_t src_idx) {
  BelowNode(dst_rotated, dst_idx) =
      BelowNode(dst_rotated, dst_idx).max(BelowNode(src_idx));
}
