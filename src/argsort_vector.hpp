// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Argsort Vectors are an associated reindexer for a reference data vector.
// This maintain the relationship between reindexer and data. Can maintain a proxy sort
// of the data vector, which can then support sorted inserts into.  Rearranging elements
// in the reindexer can avoid the heavy cost of copying heavier data objects in the
// vector. Can also be used for bidirectional or multimaps.

#pragma once

#include "reindexer.hpp"

/*
TODO: Remove this later
  AccessFunction access_fn = [](VectorType &data, const size_t i) { return data[i]; },
  ReindexFunction sort_fn =
      [](VectorType &data, const Reindexer &reindexer) {
        Reindexer::ReindexVectorInPlace<VectorType, DataType>(
            data, reindexer, data.size());
      },
*/

// Reindexer that holds a reference data vector.  The reindexer maintains a sort by
// proxy on the underlying data.

// ** Default Functions

template <typename DataType>
bool ArgsortLessThanFunction(const DataType &lhs, const DataType &rhs) {
  return lhs < rhs;
}

template <typename VectorType, typename DataType>
const DataType &ArgsortAccessFunction(VectorType data_vector, const size_t i) {
  return data_vector[i];
}

template <typename VectorType, typename DataType>
void ArgsortReindexFunction(VectorType data_vector, const Reindexer &reindexer) {
  Reindexer::ReindexVectorInPlace<VectorType, DataType>(
      data_vector, reindexer.InvertReindexer(), data_vector.size());
}

template <typename VectorType, typename DataType>
void ArgsortAddDataFunction(VectorType data_vector, VectorType data_to_add) {}

template <typename DataType, typename VectorType = std::vector<DataType>,
          typename RefVectorType = VectorType &>
class ArgsortVector {
 public:
  using AccessFunction = std::function<const DataType &(RefVectorType, const size_t)>;
  using ReindexFunction = std::function<void(RefVectorType, const Reindexer &)>;
  using LessThanFunction = std::function<bool(const DataType &, const DataType &)>;

  ArgsortVector(
      RefVectorType data_vector, std::optional<Reindexer> reindexer,
      LessThanFunction lessthan_fn = ArgsortLessThanFunction<DataType>,
      AccessFunction access_fn = ArgsortAccessFunction<RefVectorType, DataType>,
      ReindexFunction reindex_fn = ArgsortReindexFunction<RefVectorType, DataType>,
      bool is_sorted = false)
      : data_vector_(data_vector),
        is_sorted_(is_sorted),
        access_fn_(access_fn),
        reindex_fn_(reindex_fn),
        lessthan_fn_(lessthan_fn) {
    reindexer_ = reindexer.has_value()
                     ? reindexer.value()
                     : Reindexer::IdentityReindexer(data_vector.size());
    if (!is_sorted_) {
      SortReindexer();
    }
    is_sorted_ = true;
  };

  // ** Access

  bool IsSorted() const { return is_sorted_; };
  size_t Size() const { return data_vector_.size(); };

  // VectorType &GetDataVector() { return data_vector_; }
  const RefVectorType GetDataVector() const { return data_vector_; };
  // Reindexer &GetReindexer() { return reindexer_; }
  const Reindexer &GetReindexer() const { return reindexer_; };

  // Get sorted index by given unsorted index.
  size_t GetSortedIndexByUnsortedIndex(const size_t unsorted_idx) const {
    return reindexer_.GetOutputIndexByInputIndex(unsorted_idx);
  };
  // Get unsorted index by given sorted index.
  // Note: This uses linear search.
  size_t GetUnsortedIndexBySortedIndex(const size_t sorted_idx) const {
    return reindexer_.GetInputIndexByOutputIndex(sorted_idx);
  };

  // Get data by unsorted index.
  const DataType &GetDataByUnsortedIndex(const size_t unsorted_idx) const {
    return access_fn_(data_vector_, unsorted_idx);
  }
  // Get data by sorted index.
  // Note: This uses linear search.
  const DataType &GetDataBySortedIndex(const size_t sorted_idx) const {
    size_t unsorted_idx = GetUnsortedIndexBySortedIndex(sorted_idx);
    return GetDataByUnsortedIndex(unsorted_idx);
  }

  // ** Query

  template <class ForwardIt>
  ForwardIt LowerBound(ForwardIt first, ForwardIt last, const DataType &query) const {
    ForwardIt it;
    typename std::iterator_traits<ForwardIt>::difference_type count, step;
    count = std::distance(first, last);
    while (count > 0) {
      it = first;
      step = count / 2;
      std::advance(it, step);
      auto current_value = *it;
      if (lessthan_fn_(GetDataByUnsortedIndex(current_value), query)) {
        first = ++it;
        count -= step + 1;
      } else
        count = step;
    }
    return first;
  }

  template <class ForwardIt>
  ForwardIt UpperBound(ForwardIt first, ForwardIt last, const DataType &query) const {
    ForwardIt it;
    typename std::iterator_traits<ForwardIt>::difference_type count, step;
    count = std::distance(first, last);

    while (count > 0) {
      it = first;
      step = count / 2;
      std::advance(it, step);
      auto current_value = *it;
      if (!lessthan_fn_(query, GetDataByUnsortedIndex(current_value))) {
        first = ++it;
        count -= step + 1;
      } else
        count = step;
    }
    return first;
  }

  size_t FindFirstSortedIndex(const DataType query) const {
    auto lower_bound =
        LowerBound(reindexer_.GetData().begin(), reindexer_.GetData().end(), query);
    Assert(lower_bound != reindexer_.GetData().end(),
           "Searched data does not exist in data_vector.");
    return lower_bound - reindexer_.GetData().begin();
  };

  size_t FindLastSortedIndex(const DataType query) const {
    auto upper_bound =
        UpperBound(reindexer_.GetData().begin(), reindexer_.GetData().end(), query);
    Assert(upper_bound != reindexer_.GetData().end(),
           "Searched data does not exist in data_vector.");
    return upper_bound - reindexer_.GetData().begin();
  };

  // Find data in
  size_t FindUniqueSortedIndex(const DataType &data) const {
    return FindFirstSortedIndex(data);
  }

  SizePair FindRangeSortedIndex(const DataType &data) const {
    size_t range_begin = FindFirstSortedIndex(data);
    size_t range_end = FindLastSortedIndex(data);
    return {range_begin, range_end};
  };

  // Construct a sorted version of the data vector, without modifying the underlying
  // data.
  VectorType BuildSortedDataVector() const {
    VectorType sorted_vector =
        Reindexer::BuildReindexedVector(data_vector_, reindexer_);
    return sorted_vector;
  };

  // ** Modify

  // Append data_to_insert_vector to data_vector, then insert into sorted reindexer.
  void SortedInsert(RefVectorType data_to_insert_vector,
                    std::optional<bool> do_single_insert = std::nullopt) {
    if (data_to_insert_vector.empty()) {
      return;
    }
    // sort and append new data to data_vector.
    std::sort(data_to_insert_vector.begin(), data_to_insert_vector.end(), lessthan_fn_);
    data_vector_.insert(data_vector_.end(), data_to_insert_vector.begin(),
                        data_to_insert_vector.end());
    // Rough estimate -- if quantity of new data being added is more than log(N), then
    // we are better off incurring the cost of a full vector resort than doing
    // individual inserts.
    if (!do_single_insert.has_value()) {
      do_single_insert = (data_to_insert_vector.size() < log(data_vector_.size()));
    }

    if (do_single_insert.value()) {
      SizeVector sorted_output_idxs_to_add;
      for (const auto &data : data_to_insert_vector) {
        sorted_output_idxs_to_add.push_back(FindLastSortedIndex(data));
      }
      reindexer_.InsertOutputIndex(sorted_output_idxs_to_add);
    } else {
      // Add data_to_insert_vector and re-sort.
      reindexer_.AppendNextIndex(data_to_insert_vector.size());
      SortReindexer();
    }
  };

  void SortedDelete(const DataType &data_to_delete) {
    size_t id_to_delete = FindUniqueSortedIndex(data_to_delete);
    return SortedDeleteById(id_to_delete);
  };

  void SortedDelete(const RefVectorType data_to_delete_vector) {
    SizeVector ids_to_delete;
    for (size_t i = 0; i < data_to_delete_vector.size(); i++) {
      ids_to_delete.push_back(FindUniqueSortedIndex(data_to_delete_vector[i]));
    }
    return SortedDeleteById(ids_to_delete);
  };

  void SortedDeleteById(const size_t id_to_delete) {
    reindexer_.ReassignOutputIndexAndShift(id_to_delete, reindexer_.size() - 1);
  };

  void SortedDeleteById(const SizeVector &ids_to_delete) {
    for (const auto &id_to_delete : ids_to_delete) {
      SortedDeleteById(id_to_delete);
    }
  };

  // ** Transform

  // Sort reindexer using data vector ordering.
  void SortReindexer() {
    std::sort(reindexer_.GetData().begin(), reindexer_.GetData().end(),
              [this](int left, int right) -> bool {
                return lessthan_fn_(access_fn_(data_vector_, left),
                                    access_fn_(data_vector_, right));
              });
    is_sorted_ = true;
  };

  // Sort data according to the reindexer ordering.
  // Reindexer is updated to identity after sorting.
  void SortDataVector() {
    // reindex_fn_(data_vector_, reindexer_.InvertReindexer());
    reindex_fn_(data_vector_, reindexer_);
    reindexer_ = Reindexer::IdentityReindexer(data_vector_.size());
  };

  // ** Default Functions

 private:
  RefVectorType data_vector_;
  Reindexer reindexer_;
  std::optional<Reindexer> inverted_reindexer_ = std::nullopt;
  bool is_sorted_ = false;
  size_t occupancy = 0;

  AccessFunction access_fn_;
  LessThanFunction lessthan_fn_;
  ReindexFunction reindex_fn_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("ArgsortVector") {
  StringVector strings = {"d", "a", "a", "b", "g", "f", "i", "j"};
  StringVector golden_strings = StringVector(strings);
  std::sort(golden_strings.begin(), golden_strings.end());
  StringVector argsort_strings = StringVector(strings);
  ArgsortVector<std::string> argsort(argsort_strings, std::nullopt);

  // TEST_0: Check initial sort on build.
  CHECK_MESSAGE(golden_strings != strings, "TEST_0 failed.");
  CHECK_MESSAGE(strings == argsort.GetDataVector(), "TEST_0 failed.");
  CHECK_MESSAGE(golden_strings == argsort.BuildSortedDataVector(), "TEST_0 failed.");

  StringVector append_strings;
  std::string append_string;

  // TEST_1: Multiple insert in order (as).
  append_strings = {"b", "d"};
  strings.insert(strings.end(), append_strings.begin(), append_strings.end());
  golden_strings.insert(golden_strings.end(), append_strings.begin(),
                        append_strings.end());
  std::sort(golden_strings.begin(), golden_strings.end());
  argsort.SortedInsert(append_strings, true);
  CHECK_MESSAGE(golden_strings != argsort.GetDataVector(), "TEST_1 failed.");
  CHECK_MESSAGE(golden_strings == argsort.BuildSortedDataVector(), "TEST_1 failed.");

  // TEST_2: Multiple insert not in order (as batch).
  append_strings = {"a", "c", "g", "c"};
  strings.insert(strings.end(), append_strings.begin(), append_strings.end());
  golden_strings.insert(golden_strings.end(), append_strings.begin(),
                        append_strings.end());
  std::sort(golden_strings.begin(), golden_strings.end());
  argsort.SortedInsert(append_strings, false);
  CHECK_MESSAGE(golden_strings != argsort.GetDataVector(), "TEST_2 failed.");
  CHECK_MESSAGE(golden_strings == argsort.BuildSortedDataVector(), "TEST_2 failed.");

  // TEST_4: Single insert.
  append_strings = {"c"};
  strings.insert(strings.end(), append_strings.begin(), append_strings.end());
  golden_strings.insert(golden_strings.end(), append_strings.begin(),
                        append_strings.end());
  std::sort(golden_strings.begin(), golden_strings.end());
  argsort.SortedInsert(append_strings);
  CHECK_MESSAGE(golden_strings != argsort.GetDataVector(), "TEST_3 failed.");
  CHECK_MESSAGE(golden_strings == argsort.BuildSortedDataVector(), "TEST_3 failed.");

  // TEST_4: Empty insert.
  append_strings = {};
  strings.insert(strings.end(), append_strings.begin(), append_strings.end());
  golden_strings.insert(golden_strings.end(), append_strings.begin(),
                        append_strings.end());
  std::sort(golden_strings.begin(), golden_strings.end());
  argsort.SortedInsert(append_strings);
  CHECK_MESSAGE(golden_strings != argsort.GetDataVector(), "TEST_4 failed.");
  CHECK_MESSAGE(golden_strings == argsort.BuildSortedDataVector(), "TEST_4 failed.");

  // TEST_5: Single delete.
}

#endif  // DOCTEST_LIBRARY_INCLUDED
