// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Resizer is a tool for resizing padded data vectors.  After the DAG is modified, a
// resizer instance can be used to sync the resized DAG across associated data.

#pragma once

#include "sugar.hpp"

class Resizer {
 public:
  Resizer(const size_t old_count, const size_t old_spare, const size_t old_alloc,
          std::optional<const size_t> new_count, std::optional<const size_t> new_spare,
          std::optional<const size_t> explicit_alloc = std::nullopt,
          const double resizing_factor = 2.0)
      : old_count_(old_count), old_spare_(old_spare), old_alloc_(old_alloc) {
    new_count_ = (new_count.has_value() ? new_count.value() : GetOldCount());
    new_spare_ = (new_spare.has_value() ? new_spare.value() : GetOldSpare());
    new_alloc_ = GetOldAlloc();
    if (GetNewPadded() > GetOldAlloc()) {
      new_alloc_ = size_t(ceil(double(GetNewPadded()) * resizing_factor));
    }
    if (explicit_alloc.has_value()) {
      Assert(explicit_alloc.value() >= GetNewCount(),
             "Attempted to reallocate space smaller than count.");
      new_alloc_ = explicit_alloc.value() + GetOldSpare();
    }
  }

  size_t GetOldCount() const { return old_count_; }
  size_t GetNewCount() const { return new_count_; }
  size_t GetOldSpare() const { return old_spare_; }
  size_t GetNewSpare() const { return new_spare_; }
  size_t GetOldAlloc() const { return old_alloc_; }
  size_t GetNewAlloc() const { return new_alloc_; }
  size_t GetOldPadded() const { return old_count_ + old_spare_; }
  size_t GetNewPadded() const { return new_count_ + new_spare_; }

  template <typename VectorType, typename DataType>
  void ApplyResizeToEigenVector(
      VectorType &data_vec,
      std::optional<const DataType> default_val = std::nullopt) const {
    data_vec.conservativeResize(GetNewAlloc());
    data_vec.conservativeResize(GetNewPadded());
    if (default_val.has_value()) {
      for (size_t i = GetOldCount(); i < GetNewPadded(); i++) {
        data_vec[i] = default_val.value();
      }
    }
    // Fill new data with default values.
    if (default_val.has_value()) {
      for (size_t i = GetOldCount(); i < GetNewCount(); i++) {
        data_vec[i] = default_val.value();
      }
      // Fill spare data.
      for (size_t i = GetNewCount() + GetOldSpare(); i < GetNewCount() + GetNewSpare();
           i++) {
        data_vec[i] = default_val.value();
      }
    }
  }

  template <typename DataType>
  void ApplyResizeToStdVector(std::vector<DataType> &data_vec,
                              std::optional<const DataType> default_val = std::nullopt,
                              const bool copy_spare_data = false) const {
    data_vec.reserve(GetNewCount());
    data_vec.resize(GetNewPadded());
    // Move spare data from old location to new location.
    if (copy_spare_data) {
      for (size_t i = 0; i < std::min(GetOldSpare(), GetNewSpare()); i++) {
        data_vec[GetNewCount() + i] = data_vec[GetOldCount() + i];
      }
    }
    // Fill new data with default values.
    if (default_val.has_value()) {
      for (size_t i = GetOldCount(); i < GetNewCount(); i++) {
        data_vec[i] = default_val.value();
      }
      // Fill spare data.
      for (size_t i = GetNewCount() + GetOldSpare(); i < GetNewCount() + GetNewSpare();
           i++) {
        data_vec[i] = default_val.value();
      }
    }
  }

 private:
  size_t old_count_;
  size_t new_count_;
  size_t old_spare_;
  size_t new_spare_;
  size_t old_alloc_;
  size_t new_alloc_;
};
