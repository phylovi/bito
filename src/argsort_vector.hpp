// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Argsort Vectors are an associated reindexer on a reference data vector.
// The

#pragma once

#include "reindexer.hpp"

// Reindexer that holds a reference data vector.  The reindexer maintains a sort by
// proxy on the underlying data.
template <typename DataType, typename VectorType = std::vector<DataType>>
class ArgsortVector {
 public:
  ArgsortVector(
      VectorType &data_vector,
      std::function<bool(const DataType &, const DataType &)> lessthan_fn =
          [](const DataType &lhs, const DataType &rhs) { return lhs < rhs; },
      bool is_sorted = false)
      : reindexer_(Reindexer::IdentityReindexer(data_vector.size())),
        data_vector_(data_vector),
        is_sorted_(is_sorted),
        lessthan_fn_(lessthan_fn) {
    if (!is_sorted_) {
      SortReindexer();
    }
    is_sorted_ = true;
  };

  // ** Access

  bool IsSorted() const { return is_sorted_; };
  size_t Size() const { return data_vector_.size(); };

  const VectorType &GetDataVector() const { return data_vector_; };
  const Reindexer &GetReindexer() const { return reindexer_; };

  // Get sorted index by given unsorted index.
  size_t GetSortedIndexByUnsortedIndex(const size_t unsorted_idx) const {
    return reindexer_.GetOutputIndexByInputIndex(old_idx);
  };
  // Get unsorted index by given sorted index.
  // Note: This uses linear search.
  size_t GetUnsortedIndexBySortedIndex(const size_t sorted_idx) const {
    return reindexer_.GetInputIndexByOutputIndex(i);
  };

  // Get data by unsorted index.
  const DataType &GetDataByUnsortedIndex(const size_t unsorted_idx) const {
    return data_[unsorted_idx];
  }
  // Get data by sorted index.
  // Note: This uses linear search.
  const DataType &GetDataBySortedIndex(const size_t sorted_idx) const {
    size_t unsorted_idx = GetUnsortedIndexBySortedIndex(sorted_idx);
    return GetDataByUnsortedIndex(unsorted_idx);
  }

  // ** Query

  size_t FindFirstSortedIndex(const DataType &data) const {
    return std::lower_bound(
        reindexer_.begin(), reindexer_.end(), data,
        [this, &data](const size_t unsorted_idx, const DataType &data) -> bool {
          return lessthan_fn_(GetDataByUnsortedIndex(sorted_idx), )
        });
    return 0;
  };

  size_t FindLastSortedIndex(const DataType &data) const { return 0; };

  SizePair FindRangeSortedIndex(const DataType &data) const {
    size_t range_begin = FindFirstSortedIndex(data);
    size_t range_end = FindLastSortedIndex(data);
    return {range_begin, range_end};
  };

  // Construct a sorted version of the data vector, without modifying the underlying
  // data.
  VectorType BuildSortedDataVector() const {
    return Reindexer::BuildReindexedVector(data_vector_, reindexer_);
  };

  // ** Modify

  // Append data_to_insert to data_vector, and insert into sorted reindexer.
  void SortedInsert(DataType &data_to_insert) {
    // Add data element.
    data_vector_.push_back(data_to_insert);
    reindexer_.AppendNextIndex();
    // Find insert position in sorted reindexer.
    // reindexer_.GetData().insert();
    // reindexer_.GetData().insert(std::upper_bound(
    //     reindexer_.GetData().begin(), reindexer_.GetData().end(),
    //     [this](const int left, const int right) {
    //       return lessthan_fn_(reindexer_.GetData()[left],
    //       reindexer_.GetData()[right]);
    //     }));
  };

  // Append data_to_insert_vector to data_vector, then insert into sorted reindexer.
  void SortedInsert(VectorType &data_to_insert_vector) {
    // Rough estimate -- if quantity of new data being added is more than log(N), then
    // we are better off incurring the cost of a full vector resort than doing
    // individual inserts.
    if (data_to_insert_vector.size() < log(data_vector_.size())) {
      std::sort(data_to_insert_vector.begin(), data_to_insert_vector.end(),
                lessthan_fn_);
      SizeVector sorted_new_ids_to_add;
      for (const auto &data : data_to_insert_vector) {
      }
    }
    // Add data_to_insert_vector and re-sort.
    reindexer_.AppendNextIndex(data_to_insert_vector.size());
    data_vector_.insert(data_vector_.begin(), data_to_insert_vector.begin(),
                        data_to_insert_vector.end());
    SortReindexer();
  };

  void SortedDelete(DataType &data_to_delete){

  };

  void SortedDelete(VectorType &data_to_delete_vector){

  };

  void SortedDeleteById(size_t id_to_delete){};
  void SortedDeleteById(SizeVector &ids_to_delete){};

  // ** Transform

  // Sort reindexer using data vector ordering.
  void SortReindexer() {
    std::sort(reindexer_.GetData().begin(), reindexer_.GetData().end(),
              [this](int left, int right) -> bool {
                return lessthan_fn_(data_vector_[left], data_vector_[right]);
              });
    is_sorted_ = true;
  };

  // Sort data according to the reindexer ordering.
  // Reindexer is updated to identity after sorting.
  void SortDataVector() {
    Reindexer::ReindexVectorInPlace<VectorType, DataType>(data_vector_, reindexer_,
                                                          data_vector_.size());
    reindexer_ = Reindexer::IdentityReindexer(data_vector_.size());
  };

  // ** Iterator

 private:
  Reindexer reindexer_;
  std::optional<Reindexer> inverted_reindexer_ = std::nullopt;
  VectorType &data_vector_;
  bool is_sorted_ = false;
  std::function<bool(const DataType &, const DataType &)> lessthan_fn_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("ArgsortVector") {
  StringVector strings = {"d", "a", "c", "a", "b", "g", "e", "f", "i"};
  StringVector golden_strings = StringVector(strings);
  std::sort(golden_strings.begin(), golden_strings.end());
  StringVector argsort_strings = StringVector(strings);
  ArgsortVector<std::string> argsort(argsort_strings);

  CHECK_NE(golden_strings, strings);
  CHECK_EQ(strings, argsort.GetDataVector());
  CHECK_EQ(golden_strings, argsort.BuildSortedDataVector());

  StringVector append_strings = {"a", "e", "c"};
  strings.insert(strings.end(), append_strings.begin(), append_strings.end());
  golden_strings.insert(golden_strings.end(), append_strings.begin(),
                        append_strings.end());
  std::sort(golden_strings.begin(), golden_strings.end());

  argsort.SortedInsert(append_strings);

  CHECK_NE(golden_strings, argsort.GetDataVector());
  CHECK_EQ(golden_strings, argsort.BuildSortedDataVector());

  argsort.SortDataVector();

  CHECK_EQ(golden_strings, argsort.GetDataVector());
}

#endif  // DOCTEST_LIBRARY_INCLUDED
