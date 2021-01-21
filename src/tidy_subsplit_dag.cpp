// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tidy_subsplit_dag.hpp"

TidySubsplitDAG::TidySubsplitDAG() : SubsplitDAG() {}

TidySubsplitDAG::TidySubsplitDAG(EigenMatrixXb above) : above_(above){};

TidySubsplitDAG::TidySubsplitDAG(const RootedTreeCollection &tree_collection)
    : SubsplitDAG(tree_collection) {}

EigenArrayXbRef TidySubsplitDAG::BelowNode(size_t node_idx) {
  return above_.col(node_idx).array();
}

EigenArrayXb TidySubsplitDAG::AboveNode(size_t node_idx) {
  return above_.row(node_idx).array();
}
