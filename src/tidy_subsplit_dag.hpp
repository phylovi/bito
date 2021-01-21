// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A "tidy" subsplit DAG has a notion of clean and dirty vectors.

#ifndef SRC_TIDY_SUBSPLIT_DAG_HPP_
#define SRC_TIDY_SUBSPLIT_DAG_HPP_

#include "subsplit_dag.hpp"

class TidySubsplitDAG : public SubsplitDAG {
 public:
  TidySubsplitDAG();
  explicit TidySubsplitDAG(const RootedTreeCollection &tree_collection);

  EigenArrayXbRef AboveNode(size_t node_idx);
  EigenArrayXbRef BelowNode(size_t node_idx);

 private:
  // above_.(i,j) is true if i is above j, otherwise it's false.
  EigenMatrixXb above_;
};

#endif  // SRC_TIDY_SUBSPLIT_DAG_HPP_
