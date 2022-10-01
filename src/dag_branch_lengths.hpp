// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This class handles the memory and optimization for Branch Lengths.

#include "gp_dag.hpp"

template <class EigenVectorType, class DataType, class DAGElementId>
class DAGParameter {
 public:
  DAGParameter(GPDAG &dag) { dag_ = &dag; }

  // ** Counts
  size_t GetCount() { return count; }
  size_t GetSpareCount() { return spare_count; }
  size_t GetPaddedCount() { return count + spare_count_; }
  size_t GetAllocCount() { return alloc_count; }

  void SetCount(const size_t count) { count_ = count; }
  void SetSpareCount(const size_t spare_count) { spare_count_ = spare_count; }
  void SetAllocCount(const size_t alloc_count) { alloc_count_ = alloc_count; }

  // ** Access
  virtual DataType Get(const DAGElementId element_id) const {
    return data_vec_[element_id.value_];
  }
  virtual DataType GetSpare(const DAGElementId element_id) const {
    return data_vec_[element_id.value_ + GetCount()];
  }

  // ** Resize
  virtual void Resize(const size_t new_count) {
    size_t old_count = GetCount();
    SetCount(new_count);
    data_vec_.conservativeResize(GetAllocCount());
    data_vec_.conservativeResize(GetPaddedCount());
  }

 protected:
  DataType Get(const DAGElementId element_id) const {
    return data_vec_[element_id.value_];
  }

  // Data vector
  EigenVectorType data_vec_;
  // Non-ownership DAG reference.
  GPDAG *dag_ = nullptr;
  // Counts
  size_t count_ = 0;
  size_t spare_count_ = 2;
  size_t alloc_count_ = 0;
};

class DAGNodeParameter : public DAGParameter<NodeId> {
 public:
  DAGNodeParameter(GPDAG &dag) : DAGParameter<NodeId> { dag_ = &dag; }
};

class DAGEdgeParameter : public DAGParameter<EdgeId> {};

class DAGBranchLengths {
 public:
  DAGBranchLengths() {}

  DAGBranchLengths(GPDAG &dag) { dag_ = &dag; }

  EigenVectorXd &GetBranchLengths() { return branch_lengths_; }
  const EigenVectorXd &GetBranchLengths() const { return branch_lengths_; }
  const GPDAG &GetDAG() const { return *dag_; }

  void Resize() { Resize(GetDAG().EdgeCountWithLeafSubsplits()); }
  void Resize(const size_t edge_count) {}

 protected:
  //
  EigenVectorXd branch_lengths_;
  // Non-ownership DAG reference.
  GPDAG *dag_ = nullptr;
  // Counts
  size_t count_ = 0;
  size_t spare_count_ = 2;
  size_t alloc_count_ = 0;
};
