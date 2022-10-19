// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This class handles DAG General Data and DAG Branch Lengths. Handles storage and
// resizing of data vectors and branch length optimization.

#pragma once

#include "reindexer.hpp"
#include "gp_dag.hpp"

template <class VectorType, class DataType, class DAGElementId,
          size_t DataPerElement = 1>
class DAGData {
 public:
  explicit DAGData(const size_t count = 0) : data_vec_() { Resize(count); }
  DAGData(GPDAG &dag, const size_t count = 0) : data_vec_(), dag_(&dag) {
    Resize(count);
  }

  // ** Counts

  // Get data size.
  size_t GetCount() { return count_; }
  size_t GetSpareCount() { return spare_count_; }
  size_t GetPaddedCount() { return count_ + spare_count_; }
  size_t GetAllocCount() { return alloc_count_; }

  // Set data size.
  void SetCount(const size_t count) { count_ = count; }
  void SetSpareCount(const size_t spare_count) { spare_count_ = spare_count; }
  void SetAllocCount(const size_t alloc_count) { alloc_count_ = alloc_count; }

  // ** Access

  // Access element in data.
  DataType &Get(const DAGElementId element_id) { return data_vec_[element_id.value_]; }
  DataType &operator()(const DAGElementId element_id) { return Get(element_id); }
  DataType &GetSpare(const DAGElementId element_id) {
    return data_vec_[element_id.value_ + GetCount()];
  }

  // Check if there is an assigned reference DAG.
  const bool HasDAG() const { return dag_ != nullptr; }
  // Access reference DAG.
  const GPDAG &GetDAG() const { return *dag_; }

  // Get full data vector.
  VectorType &GetData() { return data_vec_; }
  const VectorType &GetData() const { return data_vec_; }

  // Set value newly added elements to be set to.
  void SetDefaultValue(const DataType default_value) { default_value_ = default_value; }

  // ** Maintanence

  virtual void Resize(std::optional<const Reindexer> reindexer = std::nullopt,
                      std::optional<const size_t> explicit_alloc = std::nullopt) {
    Failwith("Cannot call base Resize without explicit count.");
  }
  virtual void Resize(const size_t new_count,
                      std::optional<const Reindexer> reindexer = std::nullopt,
                      std::optional<const size_t> explicit_alloc = std::nullopt) {
    size_t old_count = GetCount();
    SetCount(new_count);
    if (GetPaddedCount() > GetAllocCount()) {
      SetAllocCount(size_t(ceil(double(GetPaddedCount()) * resizing_factor_)));
      if (explicit_alloc.has_value()) {
        Assert(explicit_alloc.value() >= GetCount(),
               "Attempted to reallocate space smaller than count.");
        SetAllocCount(explicit_alloc.value() + GetSpareCount());
      }
      data_vec_.conservativeResize(GetAllocCount());
    }
    data_vec_.conservativeResize(GetPaddedCount());
    if (default_value_.has_value()) {
      for (size_t i = old_count; i < GetPaddedCount(); i++) {
        data_vec_[i] = default_value_.value();
      }
    }
    if (reindexer.has_value()) {
      Reindex(reindexer.value(), new_count);
    }
  }

  virtual void Reindex(const Reindexer &reindexer, const size_t length) {
    Reindexer::ReindexInPlace<VectorType, DataType>(data_vec_, reindexer, length);
  }

 protected:
  // Data vector
  VectorType data_vec_;
  size_t data_per_element_ = DataPerElement;
  // Non-ownership DAG reference.
  GPDAG *dag_ = nullptr;
  // Value for determining growth rate of data vector.
  double resizing_factor_ = 2;
  // Value to set newly allocated data to.
  std::optional<DataType> default_value_ = std::nullopt;
  // Counts
  size_t count_ = 0;
  size_t spare_count_ = 2;
  size_t alloc_count_ = 0;
};

template <class VectorType, class DataType>
class DAGNodeData : public DAGData<VectorType, DataType, NodeId> {
 public:
  using DAGNodeDataParent = DAGData<VectorType, DataType, NodeId>;
  explicit DAGNodeData<VectorType, DataType>(const size_t count = 0)
      : DAGData<VectorType, DataType, NodeId>(count) {}
  explicit DAGNodeData<VectorType, DataType>(GPDAG &dag)
      : DAGData<VectorType, DataType, NodeId>(dag, dag.NodeCount()) {}

  virtual void Resize(const size_t new_count,
                      std::optional<const Reindexer> reindexer = std::nullopt,
                      std::optional<const size_t> explicit_alloc = std::nullopt) {
    DAGNodeDataParent::Resize(new_count, reindexer, explicit_alloc);
  }
  virtual void Resize(std::optional<const Reindexer> reindexer = std::nullopt,
                      std::optional<const size_t> explicit_alloc = std::nullopt) {
    Assert(DAGNodeDataParent::HasDAG(),
           "Must assign DAG if calling Resize without explicit count.");
    Resize(DAGNodeDataParent::GetDAG().NodeCount(), reindexer, explicit_alloc);
  }
};

template <class VectorType, class DataType>
class DAGEdgeData : public DAGData<VectorType, DataType, EdgeId> {
 public:
  using DAGEdgeDataParent = DAGData<VectorType, DataType, EdgeId>;
  explicit DAGEdgeData<VectorType, DataType>(const size_t count = 0)
      : DAGData<VectorType, DataType, EdgeId>(count) {}
  explicit DAGEdgeData<VectorType, DataType>(GPDAG &dag)
      : DAGData<VectorType, DataType, EdgeId>(dag, dag.EdgeCountWithLeafSubsplits()) {}

  virtual void Resize(const size_t new_count,
                      std::optional<const Reindexer> reindexer = std::nullopt,
                      std::optional<const size_t> explicit_alloc = std::nullopt) {
    DAGEdgeDataParent::Resize(new_count, reindexer, explicit_alloc);
  }
  virtual void Resize(std::optional<const Reindexer> reindexer = std::nullopt,
                      std::optional<const size_t> explicit_alloc = std::nullopt) {
    Assert(DAGEdgeDataParent::HasDAG(),
           "Must assign DAG if calling Resize without explicit count.");
    Resize(DAGEdgeDataParent::GetDAG().EdgeCountWithLeafSubsplits(), reindexer,
           explicit_alloc);
  }
};

using DAGNodeDoubleData = DAGNodeData<EigenVectorXd, double>;
using DAGEdgeDoubleData = DAGEdgeData<EigenVectorXd, double>;
