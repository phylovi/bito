// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Reindexers are argsorted index arrays that describe reordering of data due to
// modifications of data vectors, such as when adding NNI pairs to the SubsplitDAG.
// These reindexers can then be used to reorder associated data arrays whose indices
// correspond the ordering of SubsplitDAG data arrays, such as with the node or edge
// arrays.
//
// A reindexer is a one-to-one function that maps from an old indexing scheme to a new
// indexing scheme. In other words, if old index `i` maps to new index `j`, then
// reindexer.GetNewIndexByOldIndex(`i`) = `j`. This is implemented by an underlying
// SizeVector.  For example, if old_vector = [A, B, C] and reindexer = [1, 2, 0], then
// new_vector = [C, A, B]. Note that old_vector and reindexer must have the same size.

#pragma once

#include <numeric>

#include "eigen_sugar.hpp"
#include "sugar.hpp"

class Reindexer {
 public:
  Reindexer() : data_(){};
  Reindexer(size_t size) : data_(size){};
  Reindexer(SizeVector data) : data_(std::move(data)){};

  // ** Special Constructors
  // For each position in a identity reindexer, reindexer[`i`] = `i`.
  // E.g. for size = 5, reindexer = [0, 1, 2, 3, 4].
  static Reindexer IdentityReindexer(const size_t size);

  friend bool operator==(const Reindexer &lhs, const Reindexer &rhs) {
    return lhs.GetData() == rhs.GetData();
  }
  friend bool operator!=(const Reindexer &lhs, const Reindexer &rhs) {
    return lhs.GetData() != rhs.GetData();
  }

  size_t size() const { return data_.size(); }
  void reserve(const size_t size) { data_.reserve(size); }
  void SetReindex(const size_t old_index, const size_t new_index) {
    data_.at(old_index) = new_index;
  }
  // Find mapped new/output index corresponding to given old/input index.
  size_t GetNewIndexByOldIndex(const size_t old_index) const {
    return data_.at(old_index);
  }
  // Find mapped old/input index corresponding to new/output index. Via linear search.
  size_t GetOldIndexByNewIndex(const size_t new_index) const {
    return size_t(std::find(GetData().begin(), GetData().end(), new_index) -
                  GetData().begin());
  }
  // Add new index to end of reindexer.
  void AppendNewIndex() { data_.push_back(size()); }
  void AppendNewIndex(const size_t new_index) { data_.push_back(new_index); }
  // Get underlying index vector from reindexer.
  const SizeVector &GetData() const { return data_; }
  SizeVector &GetData() { return data_; }

  // Check if reindexer is in a valid state (contains every index exactly once,
  // ranging from 0 to reindexer_size - 1).
  bool IsValid(std::optional<size_t> length = std::nullopt) const;

  // ** Modification Operations

  // Builds new inverse reindexer of a given reindexer, such that input->output becomes
  // output->input.
  Reindexer InvertReindexer() const;

  // Builds new reindexer by removing an element identified by its index and shifting
  // other idx to maintain valid reindexer.
  Reindexer RemoveOldIndex(const size_t remove_old_idx) const;
  Reindexer RemoveNewIndex(const size_t remove_new_idx) const;

  // Builds a reindexer composing apply_reindexer onto a base_reindexer. Resulting
  // reindexer contains both reindexing operations combined.
  Reindexer ComposeWith(const Reindexer &apply_reindexer);

  // In a given reindexer, take the old_id in the reindexer and reassign it to the
  // new_id and shift over the ids strictly between old_id and new_id to ensure that the
  // reindexer remains valid. For example if old_id = 1 and new_id = 4, this method
  // would shift 1 -> 4, 4 -> 3, 3 -> 2, and 2 -> 1.
  void ReassignAndShift(const size_t old_id, const size_t new_id);

  // ** Apply Operations

  // Reindexes the given data vector according to the reindexer.
  template <typename VectorType>
  static VectorType Reindex(VectorType &old_vector, const Reindexer &reindexer,
                            std::optional<size_t> length = std::nullopt) {
    size_t reindex_size = (length.has_value() ? length.value() : old_vector.size());
    Assert(
        size_t(old_vector.size()) >= reindex_size,
        "The vector must be at least as long as reindex_size in Reindexer::Reindex.");
    Assert(size_t(reindexer.size()) >= reindex_size,
           "The reindexer must be at least as long as reindex_size in "
           "Reindexer::Reindex.");
    Assert(reindexer.IsValid(reindex_size),
           "Reindexer must be valid in Reindexer::Reindex.");
    VectorType new_vector(old_vector.size());
    // Data to reindex.
    for (size_t idx = 0; idx < reindex_size; idx++) {
      new_vector[reindexer.GetNewIndexByOldIndex(idx)] = std::move(old_vector[idx]);
    }
    // Data to copy over.
    for (size_t idx = reindex_size; idx < size_t(new_vector.size()); idx++) {
      new_vector[idx] = std::move(old_vector[idx]);
    }
    return new_vector;
  };

  // Reindexes the given data vector concatenated with additional data values.
  template <typename VectorType>
  static VectorType Reindex(VectorType &old_vector, const Reindexer &reindexer,
                            VectorType &additional_values) {
    Assert(reindexer.IsValid(), "Reindexer must be valid in Reindexer::Reindex.");
    Assert(old_vector.size() + additional_values.size() ==
               static_cast<Eigen::Index>(reindexer.size()),
           "Size of the vector and additional values must add up to the reindexer size "
           "in Reindexer::Reindex.");
    VectorType new_vector(reindexer.size());
    // Data to reindex.
    for (Eigen::Index idx = 0; idx < old_vector.size(); idx++) {
      new_vector[reindexer.GetNewIndexByOldIndex(idx)] = std::move(old_vector[idx]);
    }
    // Data to copy over.
    for (Eigen::Index idx = 0; idx < additional_values.size(); idx++) {
      new_vector[reindexer.GetNewIndexByOldIndex(old_vector.size() + idx)] =
          std::move(additional_values[idx]);
    }
    return new_vector;
  };

  // Reindex data vector in-place. Expects VectorType to have `operator[]` accessor and
  // `size()`.
  template <typename VectorType, typename DataType>
  static void ReindexInPlace(VectorType &data_vector, const Reindexer &reindexer,
                             size_t length, DataType &temp1, DataType &temp2) {
    Assert(size_t(data_vector.size()) >= length,
           "data_vector wrong size for Reindexer::ReindexInPlace.");
    Assert(size_t(reindexer.size()) >= length,
           "reindexer wrong size for Reindexer::ReindexInPlace.");
    BoolVector updated_idx = BoolVector(length, false);
    for (size_t i = 0; i < length; i++) {
      size_t old_idx = i;
      size_t new_idx = reindexer.GetNewIndexByOldIndex(i);
      if (old_idx == new_idx) {
        updated_idx[old_idx] = true;
        continue;
      }
      // Because reindexing is one-to-one function, starting at any given index in the
      // the vector, if we follow the chain of remappings from each old index to its new
      // index, we will eventually form a cycle that returns to the initial old index.
      // This avoid allocating a second data array to perform the reindex, as only two
      // temporary values are needed. Only a boolean array is needed to check for
      // already updated indexes.
      bool is_current_node_updated = updated_idx[new_idx];
      temp1 = data_vector[old_idx];
      while (is_current_node_updated == false) {
        // copy data at old_idx to new_idx, and store data at new_idx in temporary.
        temp2 = data_vector[new_idx];
        data_vector[new_idx] = temp1;
        temp1 = temp2;
        // update to next idx in cycle.
        updated_idx[new_idx] = true;
        old_idx = new_idx;
        new_idx = reindexer.GetNewIndexByOldIndex(old_idx);
        is_current_node_updated = updated_idx[new_idx];
      }
    }
  };

  // Reindex id vector. Expects VectorType to have `operator[]` accessor and
  // `size()`.
  template <typename VectorType, typename DataType>
  static void ReindexInPlace(VectorType &data_vector, const Reindexer &reindexer,
                             size_t length) {
    DataType temp1, temp2;
    Reindexer::ReindexInPlace(data_vector, reindexer, length, temp1, temp2);
  }

  // Remaps each of the ids in the vector according to the reindexer.
  template <typename DataType, typename VectorType = std::vector<DataType>>
  static void RemapIdVector(VectorType &&data_vec, const Reindexer &reindexer) {
    Assert(reindexer.IsValid(), "Reindexer must be valid in Reindexer::RemapIdVector.");
    for (const auto id : data_vec) {
      Assert(size_t(id) < reindexer.size(),
             "The vector cannot contain an id out of bounds of the reindexer in "
             "Reindexer::RemapIdVector.");
    }
    std::transform(data_vec.begin(), data_vec.end(), data_vec.begin(),
                   [reindexer](DataType &old_idx) {
                     return DataType(reindexer.GetNewIndexByOldIndex(size_t(old_idx)));
                   });
  };

  // ** I/O

  friend std::ostream &operator<<(std::ostream &os, const Reindexer &reindexer) {
    os << reindexer.GetData();
    return os;
  };

 private:
  SizeVector data_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Reindexer: IdentityReindexer") {
  // Check that IdentityReindexer returns correctly.
  Reindexer correct_default({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
  CHECK_EQ(correct_default, Reindexer::IdentityReindexer(10));
}

TEST_CASE("IsValidReindexer") {
  // Index appears more than once.
  CHECK_FALSE(Reindexer({1, 3, 0, 0}).IsValid());
  // Missing an index and/or index is out of range.
  CHECK_FALSE(Reindexer({1, 3, 4, 2}).IsValid());
  // Valid reindexer.
  CHECK(Reindexer({1, 3, 0, 2}).IsValid());
}

TEST_CASE("Reindexer: Reindex") {
  // Check that Reindex throws if the given vector and reindexer have different
  // sizes.
  SizeVector old_size_vector{7, 8, 9};
  Reindexer reindexer({2, 0, 3, 1});
  CHECK_THROWS(Reindexer::Reindex(old_size_vector, reindexer));
  // Check that Reindex returns correctly.
  reindexer = Reindexer({2, 0, 1});
  SizeVector new_size_vector = Reindexer::Reindex(old_size_vector, reindexer);
  SizeVector correct_new_size_vector{8, 9, 7};
  CHECK_EQ(new_size_vector, correct_new_size_vector);
  // Check that Reindex also works with EigenVectorXd and additional values.
  EigenVectorXd old_eigen_vector(3);
  old_eigen_vector << 7, 8, 9;
  EigenVectorXd additional_values(2);
  additional_values << 10, 11;
  reindexer = Reindexer({2, 4, 0, 3, 1});
  EigenVectorXd new_eigen_vector =
      Reindexer::Reindex(old_eigen_vector, reindexer, additional_values);
  EigenVectorXd correct_new_eigen_vector(5);
  correct_new_eigen_vector << 9, 11, 7, 10, 8;
  CHECK_EQ(new_eigen_vector, correct_new_eigen_vector);
}

TEST_CASE("Reindexer: InvertReindexer") {
  // Check that inverting a vector twice results in the original vector.
  Reindexer reindexer({1, 3, 0, 2});
  Reindexer correct_inverted_reindexer({2, 0, 3, 1});
  Reindexer inverted_reindexer = reindexer.InvertReindexer();
  CHECK_EQ(inverted_reindexer, correct_inverted_reindexer);
  Reindexer correct_reindexer = Reindexer({1, 3, 0, 2});
  reindexer = inverted_reindexer.InvertReindexer();
  CHECK_EQ(reindexer, correct_reindexer);
}

TEST_CASE("Reindexer: RemapIdVector") {
  // Check that Reindex throws if the given vector has an index out of bounds of the
  // reindexer.
  SizeVector size_vector{3, 5};
  Reindexer reindexer({2, 0, 3, 1});
  CHECK_THROWS(Reindexer::RemapIdVector<size_t>(size_vector, reindexer));
  // Check that RemapIdVector returns correctly.
  size_vector = {3, 5};
  reindexer = Reindexer({2, 0, 3, 1, 6, 4, 5});
  Reindexer::RemapIdVector<size_t>(size_vector, reindexer);
  SizeVector correct_size_vector = {1, 4};
  CHECK_EQ(size_vector, correct_size_vector);
}

TEST_CASE("Reindexer: ReassignAndShift") {
  // Check that ReassignAndShift returns correctly when old_id > new_id.
  Reindexer reindexer({0, 1, 2, 3, 4, 5, 6});
  reindexer.ReassignAndShift(4, 1);
  Reindexer correct_reindexer({0, 2, 3, 4, 1, 5, 6});
  CHECK_EQ(reindexer, correct_reindexer);
  reindexer.ReassignAndShift(5, 2);
  correct_reindexer = Reindexer({0, 3, 4, 5, 1, 2, 6});
  CHECK_EQ(reindexer, correct_reindexer);
  reindexer.ReassignAndShift(1, 3);
  correct_reindexer = Reindexer({0, 2, 4, 5, 3, 1, 6});
  CHECK_EQ(reindexer, correct_reindexer);
  // Check that ReassignAndShift returns correctly when old_id = new_id.
  reindexer = Reindexer({1, 0, 4, 6, 5, 3, 2});
  reindexer.ReassignAndShift(4, 4);
  correct_reindexer = Reindexer({1, 0, 4, 6, 5, 3, 2});
  CHECK_EQ(reindexer, correct_reindexer);
  // Check that ReassignAndShift returns correctly when old_id < new_id.
  reindexer = Reindexer({6, 0, 4, 1, 5, 3, 2});
  reindexer.ReassignAndShift(1, 5);
  correct_reindexer = Reindexer({6, 0, 3, 5, 4, 2, 1});
  CHECK_EQ(reindexer, correct_reindexer);
}

TEST_CASE("Reindexer: ComposeWith") {
  // Check that identity reindexer composed with a second reindexer results in that
  // reindexer.
  Reindexer identity_reindexer, inverted_reindexer, pairswap_reindexer,
      composed_reindexer, correct_reindexer;
  identity_reindexer = Reindexer::IdentityReindexer(6);
  inverted_reindexer = Reindexer({5, 4, 3, 2, 1, 0});
  pairswap_reindexer = Reindexer({1, 0, 3, 2, 5, 4});
  composed_reindexer = identity_reindexer;
  composed_reindexer = composed_reindexer.ComposeWith(inverted_reindexer);
  composed_reindexer = composed_reindexer.ComposeWith(pairswap_reindexer);
  correct_reindexer = Reindexer({4, 5, 2, 3, 0, 1});
  CHECK_EQ(composed_reindexer, correct_reindexer);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
