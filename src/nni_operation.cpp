// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "nni_operation.hpp"

#include "bitset.hpp"
#include "subsplit_dag.hpp"

// ** NNIOperation

int NNIOperation::Compare(const NNIOperation &nni_a, const NNIOperation &nni_b) {
  auto compare_parent = Bitset::SubsplitCompare(nni_a.parent_, nni_b.parent_);
  if (compare_parent != 0) {
    return compare_parent;
  }
  auto compare_child = Bitset::SubsplitCompare(nni_a.child_, nni_b.child_);
  return compare_child;
};

int NNIOperation::Compare(const NNIOperation &nni_b) const {
  const NNIOperation &nni_a = *this;
  return Compare(nni_a, nni_b);
}

bool operator<(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) < 0;
}
bool operator<=(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) <= 0;
}
bool operator>(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) > 0;
}
bool operator>=(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) >= 0;
}
bool operator==(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) == 0;
}
bool operator!=(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) != 0;
}

NNIOperation NNIOperation::NNIOperationFromNeighboringSubsplits(
    const Bitset parent_in, const Bitset child_in,
    const bool swap_which_child_clade_with_sister, const bool which_child_of_parent) {
  // Input: Parent(X,YZ) -> Child(Y,Z).
  Bitset X = parent_in.SubsplitGetClade(!which_child_of_parent);
  // "Y" clade can be chosen arbitrarily from (Y,Z), so "Y" is chosen based on which
  // we want to swap with "X".
  Bitset Y = child_in.SubsplitGetClade(swap_which_child_clade_with_sister);
  Bitset Z = child_in.SubsplitGetClade(!swap_which_child_clade_with_sister);
  // Output: Parent(Y,XZ) -> Child(X,Z).
  Bitset parent_out = Bitset::Subsplit(Y, X | Z);
  Bitset child_out = Bitset::Subsplit(X, Z);
  return NNIOperation(parent_out, child_out);
}

NNIOperation NNIOperation::NNIOperationFromNeighboringSubsplits(
    const Bitset parent_in, const Bitset child_in,
    const bool swap_which_child_clade_with_sister) {
  bool which_clade_of_parent =
      Bitset::SubsplitIsChildOfWhichParentClade(parent_in, child_in);
  return NNIOperationFromNeighboringSubsplits(
      parent_in, child_in, swap_which_child_clade_with_sister, which_clade_of_parent);
}

NNIOperation NNIOperation::NNIOperationFromNeighboringSubsplits(
    const bool swap_which_child_clade_with_sister) const {
  return NNIOperation::NNIOperationFromNeighboringSubsplits(
      parent_, child_, swap_which_child_clade_with_sister);
}
