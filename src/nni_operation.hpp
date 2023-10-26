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
  NNIOperation() : parent_(0), child_(0){};

  NNIOperation(Bitset parent, Bitset child) : parent_(parent), child_(child) {
    focal_clade_ = Bitset::SubsplitIsChildOfWhichParentClade(parent_, child_);
  };

  NNIOperation(const std::string &parent, const std::string child)
      : NNIOperation(Bitset(parent), Bitset(child)) {}

  enum class NNIClade { ParentFocal, ParentSister, ChildLeft, ChildRight };
  static const inline size_t NNICladeCount = 4;
  class NNICladeEnum : public EnumWrapper<NNIClade, size_t, NNICladeCount,
                                          NNIClade::ParentFocal, NNIClade::ChildRight> {
  };
  enum class NNIEdge { Parent, Sister, Focal, LeftChild, RightChild };
  static const inline size_t NNIEdgeCount = 5;
  class NNIEdgeEnum : public EnumWrapper<NNIEdge, size_t, NNIEdgeCount, NNIEdge::Parent,
                                         NNIEdge::RightChild> {};
  using NNICladeArray = NNICladeEnum::Array<NNIClade>;

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

  // Produces the neighboring NNI, resulting from a clade swap between the
  // sister clade and a child clade.
  NNIOperation GetNeighboringNNI(
      const SubsplitClade child_clade_swapped_with_sister) const;
  static NNIOperation GetNeighboringNNI(
      const Bitset parent_in, const Bitset child_in,
      const SubsplitClade child_clade_swapped_with_sister,
      const SubsplitClade focal_clade);
  // If it is not known which is focal clade, it can be inferred by this overload.
  static NNIOperation GetNeighboringNNI(
      const Bitset parent_in, const Bitset child_in,
      const SubsplitClade child_clade_swapped_with_sister);

  // ** Getters

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

  SubsplitClade WhichCladeIsFocal() const { return focal_clade_; }
  SubsplitClade WhichCladeIsSister() const {
    return SubsplitCladeEnum::Opposite(focal_clade_);
  }

  // Checks whether two NNIs are neighbors. That is whether one is the result of an NNI
  // operation on the other.
  static bool AreNNIOperationsNeighbors(const NNIOperation &nni_a,
                                        const NNIOperation &nni_b);

  // Given two neighboring NNIs returns which child clade can be swapped with the
  // sister clade in the `pre_nni` to produce the `post_nni`.
  static SubsplitClade WhichCladeSwapWithSisterToCreatePostNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni);

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
  std::string ToHashString() const;
  friend std::ostream &operator<<(std::ostream &os, const NNIOperation &nni);

  Bitset parent_;
  Bitset child_;
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
  auto nni_xy = nni_yz.GetNeighboringNNI(SubsplitClade::Left);
  CHECK_EQ(correct_nni_xy, nni_xy);
  // Swap X and Z
  auto nni_xz = nni_yz.GetNeighboringNNI(SubsplitClade::Right);
  CHECK_EQ(correct_nni_xz, nni_xz);

  // Relationship is known (child_in is the rotated clade of parent_in)
  auto nni_xy_2 = NNIOperation::GetNeighboringNNI(
      parent_in, child_in, SubsplitClade::Left, SubsplitClade::Right);
  CHECK_EQ(correct_nni_xy, nni_xy_2);
  CHECK_THROWS(NNIOperation::GetNeighboringNNI(parent_in, child_in, SubsplitClade::Left,
                                               SubsplitClade::Left));
};

TEST_CASE("NNIOperation: NNISet") {
  // Clades for NNI.
  Bitset X("100");
  Bitset Y("010");
  Bitset Z("001");
  // Initial Child and Parent.
  Bitset parent_in = Bitset::Subsplit(X, Y | Z);
  Bitset child_in = Bitset::Subsplit(Y, Z);
  NNIOperation nni_yz = NNIOperation(parent_in, child_in);
  auto nni_xy = nni_yz.GetNeighboringNNI(SubsplitClade::Left);
  auto nni_xz = nni_yz.GetNeighboringNNI(SubsplitClade::Right);
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

TEST_CASE("NNIOperation: NNI Clade Mapping") {
  // Clades for NNI.
  std::vector<Bitset> clades = {Bitset("100"), Bitset("010"), Bitset("001")};
  // Iterate over all possible assignments of {X,Y,Z} clades to {sister, left, right}.
  std::vector<std::array<size_t, 3>> assignments;
  for (size_t x = 0; x < 3; x++) {
    for (size_t y = 0; y < 3; y++) {
      if (y == x) {
        continue;
      }
      for (size_t z = 0; z < 3; z++) {
        if ((z == x) || (z == y)) {
          continue;
        }
        assignments.push_back({x, y, z});
      }
    }
  }
  // For each possible pre-NNI, check that NNI produces the correct mapping.
  for (const auto assign : assignments) {
    const Bitset X(clades[assign[0]]);
    const Bitset Y(clades[assign[1]]);
    const Bitset Z(clades[assign[2]]);
    const Bitset parent = Bitset::Subsplit(X, Y | Z);
    const Bitset child = Bitset::Subsplit(Y, Z);
    const NNIOperation pre_nni(parent, child);

    for (const auto which_clade_swap : SubsplitCladeEnum::Iterator()) {
      const auto post_nni = pre_nni.GetNeighboringNNI(which_clade_swap);
      const auto clade_map =
          NNIOperation::BuildNNICladeMapFromPreNNIToNNI(pre_nni, post_nni);
      for (const auto pre_clade_type :
           {NNIClade::ParentSister, NNIClade::ChildLeft, NNIClade::ChildRight}) {
        const auto post_clade_type = clade_map[pre_clade_type];
        CHECK_MESSAGE(
            pre_nni.GetClade(pre_clade_type) == post_nni.GetClade(post_clade_type),
            "NNI Clade Map did not produce a proper mapping.");
      }
    }
  }
}

#endif  // DOCTEST_LIBRARY_INCLUDED
