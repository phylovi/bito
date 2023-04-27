#pragma once

#include "../src/bitset.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Bitset") {
  Bitset bit_from_str = Bitset("00110100");
  Bitset bit_from_sizevec = Bitset({2, 3, 5}, 8);
  CHECK_EQ(bit_from_str, bit_from_sizevec);

  Bitset a("1100");

  CHECK_EQ(a[2], false);
  CHECK_EQ(a[1], true);

  Bitset build_up(4);
  build_up.set(1);
  build_up.set(3);
  CHECK_EQ(build_up, Bitset("0101"));

  Bitset strip_down(4, true);
  strip_down.reset(0);
  strip_down.reset(2);
  CHECK_EQ(strip_down, Bitset("0101"));

  CHECK_EQ(a.size(), 4);

  CHECK_EQ(Bitset("1100"), Bitset("1100"));
  CHECK_NE(Bitset("1100"), Bitset("0100"));

  CHECK_LT(Bitset("0100"), Bitset("0110"));
  CHECK_LT(Bitset("0100"), Bitset("0110"));
  CHECK_LT(Bitset("0010"), Bitset("0100"));
  CHECK_LE(Bitset("0010"), Bitset("0100"));
  CHECK_LE(Bitset("1100"), Bitset("1100"));

  CHECK_GT(Bitset("0110"), Bitset("0100"));
  CHECK_GT(Bitset("0110"), Bitset("0100"));
  CHECK_GT(Bitset("0100"), Bitset("0010"));
  CHECK_GE(Bitset("0100"), Bitset("0010"));
  CHECK_GE(Bitset("1100"), Bitset("1100"));

  CHECK_EQ((Bitset("1100") & Bitset("1010")), Bitset("1000"));
  CHECK_EQ((Bitset("1100") | Bitset("1010")), Bitset("1110"));
  CHECK_EQ((Bitset("1100") ^ Bitset("1010")), Bitset("0110"));
  CHECK_EQ(~Bitset("1010"), Bitset("0101"));
  CHECK_EQ(Bitset("101") + Bitset("011"), Bitset("101011"));
  CHECK_EQ(std::min(Bitset("1100"), Bitset("1010")), Bitset("1010"));

  a &= Bitset("0110");
  CHECK_EQ(a, Bitset("0100"));

  CHECK_EQ(a.All(), false);
  CHECK_EQ(Bitset(4, true).All(), true);
  CHECK_EQ(a.Any(), true);
  CHECK_EQ(Bitset(4, false).Any(), false);
  CHECK_EQ(a.None(), false);
  CHECK_EQ(Bitset(4, false).None(), true);

  a.flip();
  CHECK_EQ(a, Bitset("1011"));
  a.Minorize();
  CHECK_EQ(a, Bitset("0100"));
  a.Minorize();
  CHECK_EQ(a, Bitset("0100"));

  a.CopyFrom(Bitset("10"), 0, false);
  CHECK_EQ(a, Bitset("1000"));
  a.CopyFrom(Bitset("10"), 0, true);
  CHECK_EQ(a, Bitset("0100"));
  a.CopyFrom(Bitset("10"), 2, false);
  CHECK_EQ(a, Bitset("0110"));
  a.CopyFrom(Bitset("10"), 2, true);
  CHECK_EQ(a, Bitset("0101"));

  auto singleton = Bitset("0010");
  CHECK(singleton.IsSingleton());
  CHECK_EQ(*singleton.SingletonOption(), 2);

  CHECK_EQ(Bitset("0000").Count(), 0);
  CHECK_EQ(Bitset("0100").Count(), 1);
  CHECK_EQ(Bitset("011101").Count(), 4);

  CHECK_EQ(Bitset("1001").ToVectorOfSetBitsAsString(), "0,3");
  CHECK_EQ(Bitset("0000").ToVectorOfSetBitsAsString(), "");
}

TEST_CASE("Bitset: Clades, Subsplits, PCSPs") {
  auto p = Bitset("000111");
  // Subsplit: 000|111
  CHECK_EQ(p.SubsplitGetClade(SubsplitClade::Left), Bitset("000"));
  CHECK_EQ(p.SubsplitGetClade(SubsplitClade::Right), Bitset("111"));
  // Edge: 00|01|11
  CHECK_EQ(p.PCSPGetClade(PCSPClade::Sister), Bitset("00"));
  CHECK_EQ(p.PCSPGetClade(PCSPClade::Focal), Bitset("01"));
  CHECK_EQ(p.PCSPGetClade(PCSPClade::RightChild), Bitset("11"));

  CHECK_EQ(Bitset("11001010").SubsplitCladeUnion(), Bitset("1110"));

  CHECK_EQ(Bitset("10011100").SubsplitRotate(), Bitset("11001001"));
  CHECK_EQ(Bitset("010101").SubsplitToVectorOfSetBitsAsString(), "1|0,2");

  CHECK_EQ(Bitset("101010").SubsplitIsLeftChildOf(Bitset("111000")), true);
  // #350 commented out code
  // CHECK_EQ(Bitset::SubsplitIsChildOfWhichParentClade(Bitset("111000"),
  // Bitset("101010")), true);
  CHECK_EQ(Bitset("00100001").SubsplitIsRightChildOf(Bitset("11000011")), true);
  // CHECK_EQ(Bitset::SubsplitIsChildOfWhichParentClade(Bitset("000111"),
  // Bitset("101010")), false);
  CHECK_EQ(Bitset("010001").SubsplitIsLeftChildOf(Bitset("110001")), false);
  CHECK_EQ(Bitset("010001").SubsplitIsRightChildOf(Bitset("01000011")), false);
  // Should throw because Bitsets can't be divided into equal-sized clades.
  CHECK_THROWS(Bitset("11010").SubsplitIsLeftChildOf(Bitset("10101")));
  CHECK_THROWS(Bitset("11010").SubsplitIsRightChildOf(Bitset("10101")));

  CHECK_EQ(Bitset("101010").SubsplitIsRootsplit(), true);
  CHECK_EQ(Bitset("111000").SubsplitIsRootsplit(), false);
  CHECK_EQ(Bitset("11000001").SubsplitIsRootsplit(), false);

  CHECK_EQ(Bitset("011101").PCSPIsValid(), false);
  CHECK_EQ(Bitset("000111").PCSPIsValid(), false);
  CHECK_EQ(Bitset("100100").PCSPIsValid(), false);
  CHECK_EQ(Bitset("100011001").PCSPIsValid(), true);

  CHECK_EQ(Bitset("100011001").PCSPChildIsLeaf(), false);
  CHECK_EQ(Bitset("100011000").PCSPChildIsLeaf(), true);

  CHECK_EQ(Bitset("000111010").PCSPIsParentRootsplit(), false);
  CHECK_EQ(Bitset("000111000100").PCSPIsParentRootsplit(), false);
  CHECK_EQ(Bitset("101010000").PCSPIsParentRootsplit(), true);

  CHECK_EQ(Bitset("100011001").PCSPGetParentSubsplit(), Bitset("100011"));
  CHECK_EQ(Bitset("011100001").PCSPGetParentSubsplit(), Bitset("100011"));
  CHECK_EQ(Bitset("100011001").PCSPGetChildSubsplit(), Bitset("010001"));
  CHECK_EQ(Bitset("100001110001").PCSPGetChildSubsplit(), Bitset("01100001"));
  CHECK_EQ(Bitset("100001110001").PCSPGetChildSubsplitTaxonCounts(), SizePair({1, 2}));
  CHECK_EQ(Bitset("100000111100101").PCSPGetChildSubsplitTaxonCounts(),
           SizePair({2, 2}));

  CHECK_EQ(Bitset::Singleton(4, 2), Bitset("0010"));

  CHECK_EQ(Bitset("100010"), Bitset::Subsplit(Bitset("100"), Bitset("010")));
  CHECK_EQ(Bitset("110001"), Bitset::Subsplit(Bitset("001"), Bitset("110")));
  // Invalid clade pair.
  CHECK_THROWS(Bitset::Subsplit(Bitset("1100"), Bitset("001")));
  CHECK_THROWS(Bitset::Subsplit(Bitset("111"), Bitset("001")));

  CHECK_EQ(Bitset("000110010"), Bitset::PCSP(Bitset("110000"), Bitset("100010")));
  CHECK_EQ(Bitset("110001000"), Bitset::PCSP(Bitset("110001"), Bitset("001000")));
  // Invalid parent-child pair.
  CHECK_THROWS(Bitset::PCSP(Bitset("110001"), Bitset("010001")));
  CHECK_THROWS(Bitset::PCSP(Bitset("11000101"), Bitset("010001")));
  CHECK_THROWS(Bitset::PCSP(Bitset("110001"), Bitset("110100")));

  CHECK_EQ(Bitset::RootsplitSubsplitOfClade(Bitset("0011")), Bitset("11000011"));
  CHECK_EQ(Bitset::PCSPFromUCAToRootsplit(Bitset("11000011")), Bitset("000011110011"));

  CHECK_EQ(Bitset("010000").SubsplitIsLeaf(), true);
  CHECK_EQ(Bitset("010010").SubsplitIsLeaf(), false);
  CHECK_EQ(Bitset("111000").SubsplitIsLeaf(), false);
  CHECK_EQ(Bitset::LeafSubsplitOfNonemptyClade(Bitset("010")), Bitset("010000"));
  CHECK_EQ(Bitset::LeafSubsplitOfParentSubsplit(Bitset("100001")), Bitset("001000"));
  CHECK_THROWS(Bitset::LeafSubsplitOfParentSubsplit(Bitset("100011")));
  CHECK_EQ(Bitset::PCSPFromRightParentCladeToLeaf(Bitset("100001")),
           Bitset("100001000"));
  CHECK_THROWS(Bitset::PCSPFromRightParentCladeToLeaf(Bitset("0000110")));
  CHECK_THROWS(Bitset::PCSPFromRightParentCladeToLeaf(Bitset("100101")));

  // Restrict a bitset.
  CHECK_EQ(Bitset::Remap(Bitset("10101010101"), {0, 2, 4, 6, 8, 10}), Bitset("111111"));
  // If we apply this remap 3 times we should get back to where we started.
  SizeOptionVector rotate120{6, 7, 8, 0, 1, 2, 3, 4, 5};
  auto to_rotate = Bitset("110010100");
  CHECK_EQ(Bitset::Remap(Bitset::Remap(Bitset::Remap(to_rotate, rotate120), rotate120),
                         rotate120),
           to_rotate);
  // "Lift" a bitset.
  CHECK_EQ(Bitset::Remap(Bitset("11"), {0, std::nullopt, 1}), Bitset("101"));
}

TEST_CASE("Bitset: Subsplit Sort") {
  Bitset bitset_a = Bitset::Subsplit("01001", "00100");
  CHECK_MESSAGE(Bitset::SubsplitCompare(bitset_a, bitset_a) == 0,
                "Equality: bitset_a should be equal to itself");
  // Count of bitset_a (3) comes before count of bitset_b (4).
  Bitset bitset_b = Bitset::Subsplit("00100", "01011");
  CHECK_MESSAGE(
      Bitset::SubsplitCompare(bitset_a, bitset_b) < 0,
      "Bit Count: bitset_a should be smaller/earlier sorted value than bitset_b.");
  // Union of bitset_a ("01101") comes before union of bitset_c ("11100"), counts are
  // equal.
  Bitset bitset_c = Bitset::Subsplit("01000", "10100");
  CHECK_MESSAGE(
      Bitset::SubsplitCompare(bitset_a, bitset_c) < 0,
      "Union: bitset_a should be smaller/earlier sorted value than bitset_c.");
  // Sorted clade of bitset_a ("01001") comes before sorted clade of bitset_d ("01100"),
  // counts and unions are equal.
  Bitset bitset_d = Bitset::Subsplit("00001", "01100");
  CHECK_MESSAGE(
      Bitset::SubsplitCompare(bitset_a, bitset_d) < 0,
      "Sorted Clade: bitset_a should be smaller/earlier sorted value than bitset_d.");
}

#endif  // DOCTEST_LIBRARY_INCLUDED
