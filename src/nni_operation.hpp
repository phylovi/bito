// Copyright 2019-2022 bito project contributors.
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
// The `NNISet` is a set of `NNIOperations` used to account for all "adjacent" NNIs
// to a SubsplitDAG.  That is, all output parent/child pairs which can be generated from
// a single NNI operation on an input parent/child pair taken from the set of all the
// parent/child pairs currently in the SubsplitDAG, where the output parent/child pair
// is not also already in the SubpslitDAG.

#pragma once

#include "bitset.hpp"
#include "sugar.hpp"

// * Nearest Neighbor Interchange Operation
// NNIOperation stores output parent/child pair which are the product of an NNI.
class NNIOperation {
 public:
  NNIOperation() : parent_(0), child_(0){};
  NNIOperation(Bitset parent, Bitset child) : parent_(parent), child_(child){};

  // ** Comparator
  // NNIOperations are ordered according to the std::bitset ordering of their parent
  // subsplit, then the std::bitset order their child subsplit.
  static int Compare(const NNIOperation &nni_a, const NNIOperation &nni_b);
  int Compare(const NNIOperation &nni_b) const;

  friend bool operator<(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator>(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator<=(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator>=(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator==(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator!=(const NNIOperation &lhs, const NNIOperation &rhs);

  // ** Special Constructors:
  // Produces the output NNIOperation from an input subsplit pair that results from an
  // NNI Swap according to `swap_which_child_clade_with_sister`.
  NNIOperation NNIOperationFromNeighboringSubsplits(
      const bool swap_which_child_clade_with_sister) const;
  static NNIOperation NNIOperationFromNeighboringSubsplits(
      const Bitset parent_in, const Bitset child_in,
      const bool swap_which_child_clade_with_sister, const bool which_child_of_parent);
  // If it is not known whether the child is sorted/rotated, it can be inferred by
  // this overload.
  static NNIOperation NNIOperationFromNeighboringSubsplits(
      const Bitset parent_in, const Bitset child_in,
      const bool swap_which_child_clade_with_sister);

  // ** Miscellaneous

  bool IsValid() { return Bitset::SubsplitIsParentChildPair(parent_, child_); };

  friend std::ostream &operator<<(std::ostream &os, const NNIOperation &nni) {
    os << "{" << nni.parent_.SubsplitToString() << ", " << nni.child_.SubsplitToString()
       << "}";
    return os;
  };

  Bitset parent_;
  Bitset child_;
};

using NNISet = std::set<NNIOperation>;

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
  NNIOperation nni_yz = NNIOperation(parent_in, child_in);
  // Correct Solutions.
  Bitset correct_parent_xy = Bitset::Subsplit(Y, X | Z);
  Bitset correct_child_xy = Bitset::Subsplit(X, Z);
  NNIOperation correct_nni_xy = NNIOperation(correct_parent_xy, correct_child_xy);
  Bitset correct_parent_xz = Bitset::Subsplit(Z, Y | X);
  Bitset correct_child_xz = Bitset::Subsplit(Y, X);
  NNIOperation correct_nni_xz = NNIOperation(correct_parent_xz, correct_child_xz);

  // Swap X and Y
  auto nni_xy = nni_yz.NNIOperationFromNeighboringSubsplits(false);
  CHECK_EQ(correct_nni_xy, nni_xy);
  // Swap X and Z
  auto nni_xz = nni_yz.NNIOperationFromNeighboringSubsplits(true);
  CHECK_EQ(correct_nni_xz, nni_xz);

  // Relationship is known (child_in is the rotated clade of parent_in)
  auto nni_xy_2 = NNIOperation::NNIOperationFromNeighboringSubsplits(
      parent_in, child_in, false, true);
  CHECK_EQ(correct_nni_xy, nni_xy_2);
  CHECK_THROWS(NNIOperation::NNIOperationFromNeighboringSubsplits(parent_in, child_in,
                                                                  false, false));
};

TEST_CASE("NNISet") {
  // Clades for NNI.
  Bitset X("100");
  Bitset Y("010");
  Bitset Z("001");
  // Initial Child and Parent.
  Bitset parent_in = Bitset::Subsplit(X, Y | Z);
  Bitset child_in = Bitset::Subsplit(Y, Z);
  NNIOperation nni_yz = NNIOperation(parent_in, child_in);
  auto nni_xy = nni_yz.NNIOperationFromNeighboringSubsplits(false);
  auto nni_xz = nni_yz.NNIOperationFromNeighboringSubsplits(true);
  // Insert NNIs in various orders.
  NNISet set_of_nnis_1 = NNISet();
  set_of_nnis_1.insert(nni_yz);
  set_of_nnis_1.insert(nni_xy);
  set_of_nnis_1.insert(nni_xz);
  NNISet set_of_nnis_2 = NNISet();
  set_of_nnis_2.insert(nni_xy);
  set_of_nnis_2.insert(nni_xz);
  set_of_nnis_2.insert(nni_yz);
  // Check proper ordering.
  for (const auto &set_of_nnis : {set_of_nnis_1, set_of_nnis_2}) {
    NNIOperation prv_nni = *set_of_nnis.begin();
    for (const auto &nni : set_of_nnis) {
      CHECK_MESSAGE(nni >= prv_nni, "NNIs not ordered in NNISet.");
    }
  }
}

#endif  // DOCTEST_LIBRARY_INCLUDED
