// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Reindexers are argsorted index arrays that describe reordering of data due to
// modifications of data vectors, such as when adding NNI pairs to the SubsplitDAG.
// These reindexers can then be used to reorder associated data arrays whose indices
// correspond the ordering of SubsplitDAG data arrays, such as with the node or edge
// arrays.
//
// A reindexer is a one-to-one function that maps from an input indexing scheme to an
// output indexing scheme. In other words, if input index `i` maps to output index `j`,
// then reindexer.GetOutputIndexByInputIndex(`i`) = `j`. This is implemented by an
// underlying SizeVector.  For example, if input_vector = [A, B, C] and reindexer = [1,
// 2, 0], then output_vector = [C, A, B]. Note that input_vector and reindexer must have
// the same size.

#pragma once

#include <numeric>
#include <type_traits>

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

  // ** Comparator

  friend bool operator==(const Reindexer &lhs, const Reindexer &rhs) {
    return lhs.GetData() == rhs.GetData();
  }
  friend bool operator!=(const Reindexer &lhs, const Reindexer &rhs) {
    return lhs.GetData() != rhs.GetData();
  }

  // ** Access

  size_t size() const { return data_.size(); }
  void reserve(const size_t size) { data_.reserve(size); }
  // Set mapping from input to output index.
  void SetReindex(const size_t input_idx, const size_t output_idx) {
    data_.at(input_idx) = output_idx;
  }
  // Find mapped new/output index corresponding to given old/input index.
  size_t GetOutputIndexByInputIndex(const size_t input_idx) const {
    return data_.at(input_idx);
  }
  // Find mapped old/input index corresponding to new/output index.
  // Note: This uses a linear search. If doing many lookups, do lookup on
  // inverted reindexer instead.
  size_t GetInputIndexByOutputIndex(
      const size_t output_idx,
      std::optional<Reindexer> inverted_reindexer = std::nullopt) const;

  // Get underlying index vector from reindexer.
  const SizeVector &GetData() const { return data_; }
  SizeVector &GetData() { return data_; }

  // ** Iterator

  // SizeVector::iterator begin() { return GetData().begin(); }
  // SizeVector::iterator end() { return GetData().end(); }
  // const SizeVector::iterator begin() const { return GetData().begin(); }
  // const SizeVector::iterator end() const { return GetData().end(); }

  // ** Modify

  // Ressigns input_idx from associated input_idx->old_output_idx mapping to
  // the input_idx->new_output_idx mapping, while maintaining relative output ordering
  // of all other input_idxs.
  // - For example, for old_output_idx = 1 and new_output_idx = 4, this will shift
  // output_idxs: 1 -> 4, 4 -> 3, 3 -> 2, and 2 -> 1.
  void ReassignOutputIndexAndShift(const size_t old_output_idx,
                                   const size_t new_output_idx);

  // Append next append_count ordered indices to reindexer.
  void AppendNextIndex(const size_t append_count = 1);
  // Append given index to reindexer.
  void AppendOutputIndex(const size_t output_idx);

  // Builds new reindexer by removing an element identified by its index and shifting
  // other idx to maintain valid reindexer.
  void RemoveInputIndex(const size_t input_idx_to_remove);
  void RemoveOutputIndex(const size_t output_idx_to_remove);
  // Remove vector of indices.
  void RemoveInputIndex(SizeVector &input_idx_to_remove);
  void RemoveOutputIndex(SizeVector &output_idx_to_remove,
                         std::optional<Reindexer> inverted_reindexer = std::nullopt);

  // Insert next available input indices into specified output indices.
  void InsertInputIndex(const size_t input_idx_to_add);
  void InsertOutputIndex(const size_t output_idx_to_add);
  // Insert next available input indices into specified output indices.
  void InsertInputIndex(SizeVector &sorted_input_idxs_to_add);
  void InsertOutputIndex(SizeVector &sorted_output_idxs_to_add);

  // ** Transform

  // Make indexer into IdentityReindexer.
  void IdentityReindexer();

  // Builds new inverse reindexer of a given reindexer, such that input->output
  // becomes output->input.
  Reindexer InvertReindexer() const;

  // Builds a reindexer composing apply_reindexer onto a base_reindexer. Resulting
  // reindexer contains both reindexing operations combined.
  Reindexer ComposeWith(const Reindexer &apply_reindexer);

  // ** Apply to Data Vector
  // Applies reindexing to a supplied data vector of VectorType.
  // Note: Expects VectorType to have `operator[]` accessor and `size()`.

  // Reindexes the given data vector according to the reindexer.
  template <typename VectorType>
  static VectorType ReindexVector(VectorType &old_vector, const Reindexer &reindexer,
                                  std::optional<size_t> length = std::nullopt) {
    size_t reindex_size = (length.has_value() ? length.value() : old_vector.size());
    Assert(size_t(old_vector.size()) >= reindex_size,
           "The vector must be at least as long as reindex_size in "
           "Reindexer::ReindexVector.");
    Assert(size_t(reindexer.size()) >= reindex_size,
           "The reindexer must be at least as long as reindex_size in "
           "Reindexer::ReindexVector.");
    Assert(reindexer.IsValid(reindex_size),
           "Reindexer must be valid in Reindexer::ReindexVector.");
    VectorType new_vector(old_vector.size());
    // Data to reindex.
    for (size_t idx = 0; idx < reindex_size; idx++) {
      new_vector[reindexer.GetOutputIndexByInputIndex(idx)] =
          std::move(old_vector[idx]);
    }
    // Data to copy over.
    for (size_t idx = reindex_size; idx < size_t(new_vector.size()); idx++) {
      new_vector[idx] = std::move(old_vector[idx]);
    }
    return new_vector;
  };

  // Reindexes the given data vector concatenated with additional data values.
  template <typename VectorType>
  static VectorType ReindexVector(VectorType &old_vector, const Reindexer &reindexer,
                                  VectorType &additional_values) {
    Assert(reindexer.IsValid(), "Reindexer must be valid in Reindexer::ReindexVector.");
    Assert(old_vector.size() + additional_values.size() ==
               static_cast<Eigen::Index>(reindexer.size()),
           "Size of the vector and additional values must add up to the reindexer size "
           "in Reindexer::ReindexVector.");
    VectorType new_vector(reindexer.size());
    // Data to reindex.
    for (Eigen::Index idx = 0; idx < old_vector.size(); idx++) {
      new_vector[reindexer.GetOutputIndexByInputIndex(idx)] =
          std::move(old_vector[idx]);
    }
    // Data to copy over.
    for (Eigen::Index idx = 0; idx < additional_values.size(); idx++) {
      new_vector[reindexer.GetOutputIndexByInputIndex(old_vector.size() + idx)] =
          std::move(additional_values[idx]);
    }
    return new_vector;
  };

  // Reindex data vector in place.
  template <typename VectorType, typename DataType>
  static void ReindexVectorInPlace(VectorType &data_vector, const Reindexer &reindexer,
                                   size_t length, DataType &temp1, DataType &temp2) {
    Assert(size_t(data_vector.size()) >= length,
           "data_vector wrong size for Reindexer::ReindexVectorInPlace.");
    Assert(size_t(reindexer.size()) >= length,
           "reindexer wrong size for Reindexer::ReindexVectorInPlace.");
    BoolVector updated_idx = BoolVector(length, false);
    for (size_t i = 0; i < length; i++) {
      size_t input_idx = i;
      size_t output_idx = reindexer.GetOutputIndexByInputIndex(i);
      bool is_current_node_updated = updated_idx[output_idx];
      if ((input_idx == output_idx) || is_current_node_updated) {
        updated_idx[input_idx] = true;
        continue;
      }
      // Because reindexing is one-to-one function, starting at any given index in the
      // the vector, if we follow the chain of remappings from each old index to its
      // new index, we will eventually form a cycle that returns to the initial old
      // index. This avoid allocating a second data array to perform the reindex, as
      // only two temporary values are needed. Only a boolean array is needed to check
      // for already updated indexes.
      temp1 = std::move(data_vector[input_idx]);
      while (is_current_node_updated == false) {
        // copy data at input_idx to output_idx, and store data at output_idx in
        // temporary.
        temp2 = std::move(data_vector[output_idx]);
        data_vector[output_idx] = std::move(temp1);
        temp1 = std::move(temp2);
        // update to next idx in cycle.
        updated_idx[output_idx] = true;
        input_idx = output_idx;
        output_idx = reindexer.GetOutputIndexByInputIndex(input_idx);
        is_current_node_updated = updated_idx[output_idx];
      }
    }
  };

  // Reindex data vector in place.
  // Overload provides its own temporary data for moving data if none provided.
  template <typename VectorType, typename DataType>
  static void ReindexVectorInPlace(VectorType &data_vector, const Reindexer &reindexer,
                                   size_t length) {
    DataType temp1, temp2;
    Reindexer::ReindexVectorInPlace<VectorType, DataType>(data_vector, reindexer,
                                                          length, temp1, temp2);
  }

  // Remaps each of the ids in the vector according to the reindexer.
  template <typename VectorType>
  static void RemapIdVector(VectorType &&vector, const Reindexer &reindexer) {
    Assert(reindexer.IsValid(), "Reindexer must be valid in Reindexer::RemapIdVector.");
    for (size_t id : vector) {
      Assert(id < reindexer.size(),
             "The vector cannot contain an id out of bounds of the reindexer in "
             "Reindexer::RemapIdVector.");
    }
    std::transform(vector.begin(), vector.end(), vector.begin(),
                   [reindexer](size_t input_idx) {
                     return reindexer.GetOutputIndexByInputIndex(input_idx);
                   });
  };

  // Builds vector reindexed according to reindexer.
  // Leaves original data vector unmodified.
  template <typename VectorType>
  static VectorType BuildReindexedVector(const VectorType &old_vector,
                                         const Reindexer &reindexer,
                                         std::optional<size_t> length = std::nullopt) {
    size_t reindex_size = (length.has_value() ? length.value() : old_vector.size());
    Assert(size_t(old_vector.size()) >= reindex_size,
           "The vector must be at least as long as reindex_size in "
           "Reindexer::ReindexVector.");
    Assert(size_t(reindexer.size()) >= reindex_size,
           "The reindexer must be at least as long as reindex_size in "
           "Reindexer::ReindexVector.");
    Assert(reindexer.IsValid(reindex_size),
           "Reindexer must be valid in Reindexer::ReindexVector.");
    VectorType new_vector(old_vector.size());
    // Data to reindex.
    for (size_t idx = 0; idx < reindex_size; idx++) {
      new_vector[idx] = old_vector[reindexer.GetOutputIndexByInputIndex(idx)];
    }
    // Data to copy over.
    for (size_t idx = reindex_size; idx < size_t(new_vector.size()); idx++) {
      new_vector[idx] = old_vector[idx];
    }
    return new_vector;
  };

  // ** Miscellaneous

  friend std::ostream &operator<<(std::ostream &os, const Reindexer &reindexer);

  // Check if reindexer is in a valid state (contains every index exactly once,
  // ranging from 0 to reindexer_size - 1).
  bool IsValid(std::optional<size_t> length = std::nullopt) const;

 private:
  // Stores reindexing as a mapping of index-to-value pairs in vector.
  // (old index = data_'s position) -> (new index = data_'s value)
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
  CHECK_THROWS(Reindexer::ReindexVector(old_size_vector, reindexer));
  // Check that Reindex returns correctly.
  reindexer = Reindexer({2, 0, 1});
  SizeVector new_size_vector = Reindexer::ReindexVector(old_size_vector, reindexer);
  SizeVector correct_new_size_vector{8, 9, 7};
  CHECK_EQ(new_size_vector, correct_new_size_vector);
  // Check that Reindex also works with EigenVectorXd and additional values.
  EigenVectorXd old_eigen_vector(3);
  old_eigen_vector << 7, 8, 9;
  EigenVectorXd additional_values(2);
  additional_values << 10, 11;
  reindexer = Reindexer({2, 4, 0, 3, 1});
  EigenVectorXd new_eigen_vector =
      Reindexer::ReindexVector(old_eigen_vector, reindexer, additional_values);
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
  CHECK_THROWS(Reindexer::RemapIdVector(size_vector, reindexer));
  // Check that RemapIdVector returns correctly.
  size_vector = {3, 5};
  reindexer = Reindexer({2, 0, 3, 1, 6, 4, 5});
  Reindexer::RemapIdVector(size_vector, reindexer);
  SizeVector correct_size_vector = {1, 4};
  CHECK_EQ(size_vector, correct_size_vector);
}

TEST_CASE("Reindexer: ReassignOutputIndexAndShift") {
  // Check that ReassignOutputIndexAndShift returns correctly when old_id > new_id.
  Reindexer reindexer({0, 1, 2, 3, 4, 5, 6});
  reindexer.ReassignOutputIndexAndShift(4, 1);
  Reindexer correct_reindexer({0, 2, 3, 4, 1, 5, 6});
  CHECK_EQ(reindexer, correct_reindexer);
  reindexer.ReassignOutputIndexAndShift(5, 2);
  correct_reindexer = Reindexer({0, 3, 4, 5, 1, 2, 6});
  CHECK_EQ(reindexer, correct_reindexer);
  reindexer.ReassignOutputIndexAndShift(1, 3);
  correct_reindexer = Reindexer({0, 2, 4, 5, 3, 1, 6});
  CHECK_EQ(reindexer, correct_reindexer);
  // Check that ReassignOutputIndexAndShift returns correctly when old_id = new_id.
  reindexer = Reindexer({1, 0, 4, 6, 5, 3, 2});
  reindexer.ReassignOutputIndexAndShift(4, 4);
  correct_reindexer = Reindexer({1, 0, 4, 6, 5, 3, 2});
  CHECK_EQ(reindexer, correct_reindexer);
  // Check that ReassignOutputIndexAndShift returns correctly when old_id < new_id.
  reindexer = Reindexer({6, 0, 4, 1, 5, 3, 2});
  reindexer.ReassignOutputIndexAndShift(1, 5);
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

TEST_CASE("Reindexer: Insert/Remove") {
  Reindexer reindexer = Reindexer({2, 3, 6, 4, 1, 5, 0});
  Reindexer reindexer_test, reindexer_golden;
  reindexer_test = Reindexer(reindexer);
  reindexer_test.RemoveOutputIndex(4);
  reindexer_golden = Reindexer({2, 3, 5, 1, 4, 0});
  CHECK_EQ(reindexer_test, reindexer_golden);
  reindexer_test = Reindexer(reindexer);
  reindexer_test.RemoveInputIndex(4);
  reindexer_golden = Reindexer({1, 2, 5, 3, 4, 0});
  CHECK_EQ(reindexer_test, reindexer_golden);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
