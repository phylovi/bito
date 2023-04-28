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

#include "sugar.hpp"
#include "bitset.hpp"
#include "subsplit_dag_storage.hpp"

// * Nearest Neighbor Interchange Operation
// NNIOperation stores output parent/child pair which are the product of an NNI.
class NNIOperation {
 public:
  NNIOperation() : parent_(0), child_(0), is_focal_clade_on_right_(){};

  NNIOperation(Bitset parent, Bitset child) : parent_(parent), child_(child) {
    const SubsplitClade clade =
        Bitset::SubsplitIsChildOfWhichParentClade(parent_, child_);
    is_focal_clade_on_right_ = (clade == SubsplitClade::Right);
  };

  NNIOperation(const std::string &parent, const std::string child)
      : NNIOperation(Bitset(parent), Bitset(child)) {}

  static const inline size_t NNICladeCount = 4;
  enum class NNIClade : size_t { ParentFocal, ParentSister, ChildLeft, ChildRight };
  using NNICladeEnum = EnumWrapper<NNIClade, size_t, NNICladeCount,
                                   NNIClade::ParentFocal, NNIClade::ChildRight>;
  using NNICladeArray = EnumArray<NNIClade, NNICladeCount, NNIClade>;

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

  // ** Special Constructors

  // Produces the output NNIOperation from an input subsplit pair that results from an
  // NNI Swap according to `swap_which_child_clade_with_sister`.
  NNIOperation NNIOperationFromNeighboringSubsplits(
      const bool is_sister_swapped_with_right_child) const;
  static NNIOperation NNIOperationFromNeighboringSubsplits(
      const Bitset parent_in, const Bitset child_in,
      const bool is_sister_swapped_with_right_child, const bool which_child_of_parent);
  // If it is not known whether the child is sorted/rotated, it can be inferred by
  // this overload.
  static NNIOperation NNIOperationFromNeighboringSubsplits(
      const Bitset parent_in, const Bitset child_in,
      const bool is_sister_swapped_with_right_child);

  // ** Getters

  bool WhichParentCladeIsFocalClade() const { return is_focal_clade_on_right_; };
  Bitset GetClade(const NNIClade &nni_clade) const {
    switch (nni_clade) {
      case NNIClade::ParentFocal:
        return GetFocalClade();
      case NNIClade::ParentSister:
        return GetSisterClade();
      case NNIClade::ChildLeft:
        return GetLeftChildClade();
      case NNIClade::ChildRight:
        return GetRightChildClade();
      default:
        Failwith("Invalid NNIClade given.");
    }
  };
  const Bitset &GetParent() const { return parent_; };
  const Bitset &GetChild() const { return child_; };
  Bitset GetCentralEdgePCSP() const { return Bitset::PCSP(GetParent(), GetChild()); }
  Bitset GetFocalClade() const {
    return GetParent().SubsplitGetClade(WhichCladeIsFocal());
  };
  Bitset GetSisterClade() const {
    return GetParent().SubsplitGetClade(WhichCladeIsSister());
  }
  Bitset GetLeftChildClade() const {
    return GetChild().SubsplitGetClade(SubsplitClade::Left);
  }
  Bitset GetRightChildClade() const {
    return GetChild().SubsplitGetClade(SubsplitClade::Right);
  }

  // ** Query

  SubsplitClade WhichCladeIsFocal() const {
    return (is_focal_clade_on_right_ ? SubsplitClade::Right : SubsplitClade::Left);
  }
  SubsplitClade WhichCladeIsSister() const {
    return (is_focal_clade_on_right_ ? SubsplitClade::Left : SubsplitClade::Right);
  }

  // Checks whether two NNIs are neighbors. That is whether one is the result of an NNI
  // operation on the other.
  static bool AreNNIOperationsNeighbors(const NNIOperation &nni_a,
                                        const NNIOperation &nni_b);

  // #350: use SubsplitClade instead of bools
  // Given two neighboring NNIs returns which child clade can be swapped with the
  // sister clade in the `pre_nni` to produce the `post_nni`.
  static bool WhichCladeSwapWithSisterToCreatePostNNI(const NNIOperation &pre_nni,
                                                      const NNIOperation &post_nni);

  // ** Miscellaneous

  size_t Hash() const { return GetParent().Hash() & GetChild().Hash(); }

  // Finds mappings of sister, left child, and right child clades from Pre-NNI to NNI.
  static NNICladeArray BuildNNICladeMapFromPreNNIToNNI(const NNIOperation &pre_nni,
                                                       const NNIOperation &post_nni);

  std::pair<Direction, SubsplitClade> GetDirectionAndSubsplitCladeByNNIClade(
      const NNIClade &nni_clade) {
    switch (nni_clade) {
      case NNIClade::ParentFocal:
        return {Direction::Rootward, WhichCladeIsFocal()};
      case NNIClade::ParentSister:
        return {Direction::Rootward, WhichCladeIsSister()};
      case NNIClade::ChildLeft:
        return {Direction::Leafward, SubsplitClade::Left};
      case NNIClade::ChildRight:
        return {Direction::Leafward, SubsplitClade::Right};
      default:
        Failwith("Invalid NNIClade given.");
    }
  }
  // Checks that NNI Operation is in valid state.
  // - Parent and Child are adjacent Subsplits.
  bool IsValid();

  std::string ToString() const;
  friend std::ostream &operator<<(std::ostream &os, const NNIOperation &nni);

  Bitset parent_;
  Bitset child_;
  bool is_focal_clade_on_right_;
  SubsplitClade focal_clade_;
};

using NNISet = std::set<NNIOperation>;
using NNIVector = std::vector<NNIOperation>;
using NNIClade = NNIOperation::NNIClade;
using NNICladeEnum = NNIOperation::NNICladeEnum;

// This is how we inject a hash routine and a custom comparator into the std
// namespace so that we can use unordered_map and unordered_set.
// https://en.cppreference.com/w/cpp/container/unordered_map
namespace std {
template <>
struct hash<NNIOperation> {
  size_t operator()(const NNIOperation &x) const { return x.Hash(); }
};
template <>
struct equal_to<NNIOperation> {
  bool operator()(const NNIOperation &lhs, const NNIOperation &rhs) const {
    return lhs == rhs;
  }
};
}  // namespace std
