// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Subsplit DAG NNI (Nearest Neighbor Interchange):
//
// An `NNIOperation` contains an output parent/child Subsplit pair which is the
// result of an NNI operation on an input parent/child pair. An NNI operation can be
// seen as a swapping of the branches in the SubsplitDAG, or alternatively as a a
// reordering of the set of clades in an input parent/child pair: the parent's sister
// clade, the child's left clade, and the child's right clade. For any given
// parent/child pair, there are two possible NNIs: swapping the sister clade with the
// left child clade, or swapping the sister clade with the right child clade.
//
// The `SetOfNNIs` is a set of `NNIOperations` used to account for all "adjacent" NNIs
// to a SubsplitDAG.  That is, all output parent/child pairs which can be generated from
// a single NNI operation on an input parent/child pair taken from the set of all the
// parent/child pairs currently in the SubsplitDAG, where the output parent/child pair
// is not also already in the SubpslitDAG.

#pragma once

#include "bitset.hpp"
#include "sugar.hpp"

// Forward declaration for SetOfNNIs methods.
class SubsplitDAG;

// Nearest Neighbor Interchange Operations
// NNIOperation stores output parent/child pair which are the product of an NNI.
class NNIOperation {
 public:
  NNIOperation(Bitset parent, Bitset child) : parent_(parent), child_(child){};

  // Comparator:
  // NNIOperations are ordered according to the std::bitset ordering of their parent
  // subsplit, then the std::bitset order their child subsplit.
  static int Compare(const NNIOperation &nni_a, const NNIOperation &nni_b);
  int Compare(const NNIOperation &nni_b) const;

  friend bool operator<(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator>(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator==(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator!=(const NNIOperation &lhs, const NNIOperation &rhs);

  // Special Constructors:
  // Produces the output NNIOperation from an input subsplit pair that results from an
  // NNI Swap according to `swap_which_child_clade_with_sister`.
  static NNIOperation NNIOperationFromNeighboringSubsplits(
      const Bitset parent_in, const Bitset child_in,
      const bool swap_which_child_clade_with_sister, const bool which_child_of_parent);
  // If it is not known whether the child is sorted/rotated, it can be inferred by
  // this overload.
  static NNIOperation NNIOperationFromNeighboringSubsplits(
      const Bitset parent_in, const Bitset child_in,
      const bool swap_which_child_clade_with_sister);

  Bitset parent_;
  Bitset child_;
};

// SetOfNNIs contain all NNI output parent/child pairs which are "adjacent" to the
// current SubsplitDAG. That is, pairs which are the result of an NNI on an input bitset
// pair that are currently present in the SubsplitDAG.
class SetOfNNIs {
 public:
  SetOfNNIs() : set_(){};

  friend bool operator==(const SetOfNNIs &lhs, const SetOfNNIs &rhs);
  friend bool operator!=(const SetOfNNIs &lhs, const SetOfNNIs &rhs);

  void Insert(NNIOperation nni_op);
  void Insert(Bitset parent, Bitset child);
  void Erase(NNIOperation nni_op);
  void Erase(Bitset parent, Bitset child);
  void Clear();
  size_t GetSize() const;

 private:
  std::set<NNIOperation> set_;
};

// ** NNIEvaluationEngine Methods - Issue #372
//
// Maintainence Methods: These maintain SetOfNNIs to stay consistent with the state of
// associated DAG.
//
// Freshly synchonizes SetOfNNIs to match the current state of its DAG. Wipes old NNI
// data and finds all all parent/child pairs adjacent to DAG by iterating over all
// internal edges in DAG.
void SyncSetOfNNIsWithDAG(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag);
// Updates NNI Set after given parent/child node pair have been added to the DAG.
// Removes pair from NNI Set and adds adjacent pairs coming from newly created edges.
void UpdateSetOfNNIsAfterDAGAddNodePair(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag,
                                        const Bitset &parent_bitset,
                                        const Bitset &child_bitset);
// Maintainence Helper Methods:
//
// Adds all NNIs from all (node_id, other_id) pairs, where other_id's are elements of
// the node_id_vector.
void AddAllNNIsFromNodeVectorToSetOfNNIs(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag,
                                         const size_t &node_id,
                                         const SizeVector &adjacent_node_ids,
                                         const bool is_edge_rotated,
                                         const bool is_edge_leafward);
// Based on given input NNIOperation, produces the two possible output NNIOperations
// and adds those results to the NNI Set (if results are not a member of the DAG).
void SafeAddOutputNNIsToSetOfNNIs(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag,
                                  const Bitset &parent_bitset,
                                  const Bitset &child_bitset, const bool rotated);

#ifdef DOCTEST_LIBRARY_INCLUDED

// See tree diagram at:
// https://user-images.githubusercontent.com/31897211/136849710-de0dcbe3-dc2b-42b7-b3de-dd9b1a60aaf4.gif
TEST_CASE("NNIOperation") {
  // Clades for NNI.
  Bitset X("100");
  Bitset Y("010");
  Bitset Z("001");
  // Initial Child and Parent.
  Bitset parent_in = Bitset::Subsplit(X, Y | Z);
  Bitset child_in = Bitset::Subsplit(Y, Z);
  // Correct Solutions.
  Bitset correct_parent_xy = Bitset::Subsplit(Y, X | Z);
  Bitset correct_child_xy = Bitset::Subsplit(X, Z);
  NNIOperation correct_nni_xy = NNIOperation(correct_parent_xy, correct_child_xy);
  Bitset correct_parent_xz = Bitset::Subsplit(Z, Y | X);
  Bitset correct_child_xz = Bitset::Subsplit(Y, X);
  NNIOperation correct_nni_xz = NNIOperation(correct_parent_xz, correct_child_xz);

  // Swap X and Y
  auto nni_xy =
      NNIOperation::NNIOperationFromNeighboringSubsplits(parent_in, child_in, 0);
  CHECK_EQ(correct_nni_xy, nni_xy);
  // Swap X and Z
  auto nni_xz =
      NNIOperation::NNIOperationFromNeighboringSubsplits(parent_in, child_in, 1);
  CHECK_EQ(correct_nni_xz, nni_xz);

  // Relationship is known (child_in is the rotated clade of parent_in)
  auto nni_xy_2 =
      NNIOperation::NNIOperationFromNeighboringSubsplits(parent_in, child_in, 0, 1);
  CHECK_EQ(correct_nni_xy, nni_xy_2);
  CHECK_THROWS(
      NNIOperation::NNIOperationFromNeighboringSubsplits(parent_in, child_in, 0, 0));
}

#endif  // DOCTEST_LIBRARY_INCLUDED

