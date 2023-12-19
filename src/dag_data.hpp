// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This class handles DAG General Data stored on nodes or edges. Handles storage and
// resizing of data vectors.

#pragma once

#include "resizer.hpp"
#include "reindexer.hpp"
#include "gp_dag.hpp"
#include "graft_dag.hpp"

class GPDAG;
class GraftDAG;

template <class VectorType, class DataType, class DAGElementId,
          size_t DataPerElement = 1>
class DAGData {
 public:
  explicit DAGData(const std::optional<const size_t> count = std::nullopt,
                   std::optional<DataType> default_value = std::nullopt)
      : data_vec_() {
    size_t init_count = (count.has_value() ? count.value() : 0);
    size_t init_alloc = std::max(init_alloc_, (init_count + spare_count_) * 2);
    Resize(init_count, spare_count_, init_alloc);
    if (default_value.has_value()) {
      SetDefaultValue(default_value.value());
      FillWithDefault();
    }
  }

  // ** Counts

  // Set data size.
  void SetCount(const size_t count) { count_ = count; }
  void SetSpareCount(const size_t spare_count) { spare_count_ = spare_count; }
  void SetAllocCount(const size_t alloc_count) { alloc_count_ = alloc_count; }

  // Get data size.
  size_t size() const { return data_vec_.size(); }
  size_t GetCount() const { return count_; }
  size_t GetSpareCount() const { return spare_count_; }
  size_t GetPaddedCount() const { return count_ + spare_count_; }
  size_t GetAllocCount() const { return alloc_count_; }

  // ** Access

  // Access element in data.
  DataType &Get(const DAGElementId element_id) {
    Assert(element_id.value_ <= GetPaddedCount(),
           "Attempted to access element out of range.");
    return data_vec_[element_id.value_];
  }
  const DataType &Get(const DAGElementId element_id) const {
    Assert(element_id.value_ <= GetPaddedCount(),
           "Attempted to access element out of range.");
    return data_vec_[element_id.value_];
  }
  DataType &operator()(const DAGElementId element_id) { return Get(element_id); }
  const DataType &operator()(const DAGElementId element_id) const {
    return Get(element_id);
  }
  DataType &GetSpare(const DAGElementId element_id) {
    return Get(DAGElementId(element_id.value_ + GetCount()));
  }
  const DataType &GetSpare(const DAGElementId element_id) const {
    return Get(DAGElementId(element_id.value_ + GetCount()));
  }

  // Get full data vector.
  VectorType &GetData() { return data_vec_; }
  const VectorType &GetData() const { return data_vec_; }
  // Get data corresponding to the DAG elements.
  VectorType GetDAGData() const { return data_vec_.segment(0, GetCount()); }
  // Get data corresponding to the DAG elements, including spare elements.
  VectorType GetPaddedDAGData() const { return data_vec_.segment(0, GetPaddedCount()); }
  // Get data corresponding to the DAG elements.
  void SetDAGData(const VectorType &data_vec) {
    Assert(GetCount() == size_t(data_vec.size()),
           "Data Vector must be of equal size to current count.");
    for (size_t i = 0; i < GetCount(); i++) {
      data_vec_[i] = data_vec[i];
    }
  }
  // Set data corresponding to the DAG elements, including spare elements.
  void SetPaddedDAGData(const VectorType &data_vec) {
    Assert(GetPaddedCount() == size_t(data_vec.size()),
           "Data Vector must be of equal size to current count.");
    for (size_t i = 0; i < GetCount(); i++) {
      data_vec_[i] = data_vec[i];
    }
  }

  // Value to assign to newly added elements.
  std::optional<DataType> HasDefaultValue() const { return default_value_.has_value(); }
  DataType GetDefaultValue() const {
    Assert(default_value_.has_value(), "Cannot be called when DefaultValue not set.");
    return default_value_.value();
  }
  void SetDefaultValue(const DataType default_value) { default_value_ = default_value; }

  // Fill all cells with given value.
  void Fill(const DataType fill_value) {
    for (size_t i = 0; i < GetPaddedCount(); i++) {
      data_vec_[i] = fill_value;
    }
  }
  void FillWithDefault() {
    Assert(default_value_.has_value(), "Cannot be called when DefaultValue not set.");
    Fill(default_value_.value());
  }

  // ** Maintanence

  // Resizes data vector.
  void Resize(std::optional<const size_t> new_count = std::nullopt,
              std::optional<const size_t> new_spare = std::nullopt,
              std::optional<const size_t> explicit_alloc = std::nullopt,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    const Resizer resizer(GetCount(), GetSpareCount(), GetAllocCount(), new_count,
                          new_spare, explicit_alloc, resizing_factor_);
    Resize(resizer, reindexer);
  }

  // Resizes data vector.
  void Resize(const Resizer &resizer,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    resizer.ApplyResizeToEigenVector<VectorType, DataType>(data_vec_, default_value_);
    SetCount(resizer.GetNewCount());
    SetSpareCount(resizer.GetNewSpare());
    SetAllocCount(resizer.GetNewAlloc());
    if (reindexer.has_value()) {
      Reindex(reindexer.value(), resizer.GetNewCount());
    }
  }

  // Reindex data vector. If reindexing during a resize, then length should be the old
  // count.
  void Reindex(const Reindexer &reindexer, const size_t length) {
    Reindexer::ReindexInPlace<VectorType, DataType>(data_vec_, reindexer, length);
  }

  auto begin() { return data_vec_.begin(); }
  auto end() { return data_vec_.end(); }

 protected:
  // Data vector
  VectorType data_vec_;
  size_t data_per_element_ = DataPerElement;
  // Value for determining growth rate of data vector.
  double resizing_factor_ = 2.0;
  // Value to set newly allocated data to.
  std::optional<DataType> default_value_ = std::nullopt;
  // Counts
  size_t count_ = 0;
  size_t spare_count_ = 10;
  size_t alloc_count_ = 0;

  // Unowned DAG reference.
  GPDAG *dag_ = nullptr;
  // Unowned GraftDAG reference.
  GraftDAG *graft_dag_ = nullptr;

  // Minimum default beginning allocated size of data vectors.
  static constexpr size_t init_alloc_ = 32;
};

template <class VectorType, class DataType>
class DAGNodeData : public DAGData<VectorType, DataType, NodeId> {
 public:
  explicit DAGNodeData<VectorType, DataType>(
      const std::optional<const size_t> count = std::nullopt,
      std::optional<DataType> default_value = std::nullopt)
      : DAGData<VectorType, DataType, NodeId>(count, default_value) {}
  explicit DAGNodeData<VectorType, DataType>(
      GPDAG &dag, std::optional<DataType> default_value = std::nullopt)
      : DAGData<VectorType, DataType, NodeId>(std::nullopt, default_value) {
    Resize(dag);
  }
  explicit DAGNodeData<VectorType, DataType>(
      GraftDAG &dag, std::optional<DataType> default_value = std::nullopt)
      : DAGData<VectorType, DataType, NodeId>(std::nullopt, default_value) {
    Resize(dag);
  }

  void Resize(const GPDAG &dag,
              std::optional<const size_t> explicit_alloc = std::nullopt,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    Resize(dag.NodeCount(), std::nullopt, explicit_alloc, reindexer);
  }

  void Resize(const GraftDAG &dag,
              std::optional<const size_t> explicit_alloc = std::nullopt,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    Resize(dag.HostNodeCount(), dag.GraftNodeCount(), explicit_alloc, reindexer);
  }

  void Resize(std::optional<const size_t> new_count = std::nullopt,
              std::optional<const size_t> new_spare = std::nullopt,
              std::optional<const size_t> explicit_alloc = std::nullopt,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    DAGData<VectorType, DataType, NodeId>::Resize(new_count, new_spare, explicit_alloc,
                                                  reindexer);
  }
};

template <class VectorType, class DataType>
class DAGEdgeData : public DAGData<VectorType, DataType, EdgeId> {
 public:
  explicit DAGEdgeData<VectorType, DataType>(
      const std::optional<const size_t> count = std::nullopt,
      std::optional<DataType> default_value = std::nullopt)
      : DAGData<VectorType, DataType, EdgeId>(count, default_value) {}
  explicit DAGEdgeData<VectorType, DataType>(
      GPDAG &dag, std::optional<DataType> default_value = std::nullopt)
      : DAGData<VectorType, DataType, EdgeId>(std::nullopt, default_value) {
    Resize(dag);
  }
  explicit DAGEdgeData<VectorType, DataType>(
      GraftDAG &dag, std::optional<DataType> default_value = std::nullopt)
      : DAGData<VectorType, DataType, EdgeId>(std::nullopt, default_value) {
    Resize(dag);
  }

  void Resize(const GPDAG &dag,
              std::optional<const size_t> explicit_alloc = std::nullopt,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    Resize(dag.EdgeCountWithLeafSubsplits(), std::nullopt, explicit_alloc, reindexer);
  }

  void Resize(const GraftDAG &dag,
              std::optional<const size_t> explicit_alloc = std::nullopt,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    Resize(dag.HostEdgeCount(), dag.GraftEdgeCount(), explicit_alloc, reindexer);
  }

  void Resize(std::optional<const size_t> new_count = std::nullopt,
              std::optional<const size_t> new_spare = std::nullopt,
              std::optional<const size_t> explicit_alloc = std::nullopt,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    DAGData<VectorType, DataType, EdgeId>::Resize(new_count, new_spare, explicit_alloc,
                                                  reindexer);
  }
};

using DAGNodeDoubleData = DAGNodeData<EigenVectorXd, double>;
using DAGEdgeDoubleData = DAGEdgeData<EigenVectorXd, double>;
using DAGNodeIntData = DAGNodeData<EigenVectorXi, int>;
using DAGEdgeIntData = DAGEdgeData<EigenVectorXi, int>;
