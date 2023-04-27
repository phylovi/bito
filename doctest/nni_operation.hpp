#pragma once

#include "../src/nni_operation.hpp"

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

TEST_CASE("NNIOperation: NNISet") {
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

    for (const auto which_clade_swap : {true, false}) {
      const auto post_nni =
          pre_nni.NNIOperationFromNeighboringSubsplits(which_clade_swap);
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
