// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tidy_subsplit_dag.hpp"

TidySubsplitDAG::TidySubsplitDAG() : SubsplitDAG() {}

TidySubsplitDAG::TidySubsplitDAG(EigenMatrixXb above) : above_(above){};

TidySubsplitDAG::TidySubsplitDAG(const RootedTreeCollection &tree_collection)
    : SubsplitDAG(tree_collection) {}

EigenArrayXbRef TidySubsplitDAG::AboveNode(size_t node_idx) {
  return above_.col(node_idx).array();
}

EigenArrayXbRef TidySubsplitDAG::BelowNode(size_t node_idx) {
  // This block is a way of getting the node_idx'th row.
  return above_.block(node_idx, 0, 1, above_.cols()).array();
  // I don't know why we can't do this:
  // return above_.row(node_idx).array();
}
