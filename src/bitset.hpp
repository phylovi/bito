// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A class to store bitsets in their various forms.
//
// This file started life as the RbBitSet class from RevBayes by Sebastian
// Hoehna, but has evolved a lot since then!
//
// For basic operations, I'm trying to follow the interface of std::bitset, though
// this class goes way beyond what std::bitset offers.
// Note that we can't use std::bitset because we don't know the size of the
// bitsets at compile time.

#ifndef SRC_BITSET_HPP_
#define SRC_BITSET_HPP_

#include <algorithm>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "sugar.hpp"

class Bitset {
 public:
  explicit Bitset(std::vector<bool> value);
  explicit Bitset(size_t n, bool initial_value = false);
  explicit Bitset(std::string);

  bool operator[](size_t i) const;
  size_t size() const;

  void set(size_t i, bool value = true);
  void reset(size_t i);
  void flip();

  bool operator==(const Bitset &x) const;
  bool operator!=(const Bitset &x) const;
  bool operator<(const Bitset &x) const;
  bool operator<=(const Bitset &x) const;
  bool operator>(const Bitset &x) const;
  bool operator>=(const Bitset &x) const;

  Bitset operator&(const Bitset &x) const;
  Bitset operator|(const Bitset &x) const;
  Bitset operator^(const Bitset &x) const;
  Bitset operator~() const;
  Bitset operator+(const Bitset &x) const;

  void operator&=(const Bitset &other);
  void operator|=(const Bitset &other);

  // These methods aren't in the bitset interface, so they get our usual PascalCase
  // naming convention.
  void Zero();
  size_t Hash() const;
  std::string ToString() const;
  // Are all of the bits 1?
  bool All() const;
  // Are any of the bits 1?
  bool Any() const;
  // Is exactly one of the bits 1?
  bool IsSingleton() const;
  // Take the minimum of the bitset and its complement.
  void Minorize();
  // Copy all of the bits from another bitset into this bitset, starting at
  // begin, and optionally flipping the bits as they get copied.
  void CopyFrom(const Bitset &other, size_t begin, bool flip);
  // If the bitset only has one bit on, then we return the location of that bit.
  // Otherwise, return nullopt.
  std::optional<uint32_t> SingletonOption() const;
  // Get the number of bits set to be 1;
  size_t Count() const;
  std::string ToIndexSetString() const;

  // These methods require the bitset to be a "subsplit bitset" of even length,
  // consisting of two equal sized "chunks" representing the two sides of the
  // subsplit.
  //
  // Flip the order of the two sides of a subsplit.
  Bitset RotateSubsplit() const;
  size_t SubsplitChunkSize() const;
  // Get the ith chunk of the subsplit.
  Bitset SubsplitChunk(size_t i) const;
  // Return a string of the PCSP in the specified number of chunks, with each chunk
  // separated by a "|".
  std::string ToStringChunked(size_t chunk_count) const;
  std::string SubsplitToString() const;
  std::string SubsplitToIndexSetString() const;
  // Is the right hand chunk all zero?
  bool SubsplitIsFake() const;

  // These functions require the bitset to be a "PCSP bitset" with three
  // equal-sized "chunks".
  // The first chunk represents the sister clade, the second the focal clade,
  // and the third chunk describes the children. The children are well defined
  // relative to the cut parent: the other part of the subsplit is the cut
  // parent setminus the child.
  // For example, `100011001` is composed of the chunks `100`, `011` and `001`.
  // If the taxa are x0, x1, and x2 then this means the parent subsplit is (A,
  // BC), and the child subsplit is (B,C). See the unit tests at the bottom for
  // more examples.
  std::string PCSPToString() const;
  bool PCSPIsValid() const;
  bool PCSPIsFake() const;
  // Do the sister and focal clades union to the whole taxon set?
  bool PCSPIsRootsplit() const;
  size_t PCSPChunkSize() const;
  Bitset PCSPChunk(size_t i) const;
  // Get the first 2/3rds of the PCSP.
  Bitset PCSPParent() const;
  // Get the second 2/3rds of the PCSP.
  Bitset PCSPWithoutParent() const;
  // Get the representation of the child subsplit in the simple form of a pair
  // of membership indicators, concatenated together into one bitset.
  Bitset PCSPChildSubsplit() const;
  // Get the number of taxa in each side of the child subsplit.
  SizePair PCSPChildSubsplitTaxonCounts() const;

  // ** Static methods
  // Make a bitset with only the specified entry turned on.
  static Bitset Singleton(size_t n, size_t which_on);
  // Make the full subsplit out of the parent subsplit, whose second half is
  // split by child_half.
  static Bitset ChildSubsplit(const Bitset &parent_subsplit, const Bitset &child_half);
  // Build a PCSP bitset out of a compatible parent-child pair of bitsets.
  static Bitset PCSPOfPair(const Bitset &parent_subsplit, const Bitset &child_subsplit,
                           bool assert_validity = true);
  // Make a "fake" subsplit, which pads the nonzero contents on the right with zero to
  // have double the width.
  static Bitset FakeSubsplit(const Bitset &nonzero_contents);
  // Make a "fake" child subsplit of a given parent; assert that the left-hand chunk of
  // the parent subsplit is non-empty and that the right-hand chunk is a singleton.
  // This fake subsplit has this singleton on the left and all zeroes on the right.
  static Bitset FakeChildSubsplit(const Bitset &parent_subsplit);
  // Make a "fake" PCSP of a given parent subsplit; assert that the left-hand chunk of
  // the parent subsplit is non-empty and that the right-hand chunk is a singleton.
  // This fake subsplit has parent subsplit on the left and all zeroes on the right.
  static Bitset FakePCSP(const Bitset &parent_subsplit);

 private:
  std::vector<bool> value_;

  static void AssertIsDisjointUnion(const Bitset &should_be_union, const Bitset &set_1,
                                    const Bitset &set_2);
};

// This is how we inject a hash routine and a custom comparator into the std
// namespace so that we can use unordered_map and unordered_set.
// https://en.cppreference.com/w/cpp/container/unordered_map
namespace std {
template <>
struct hash<Bitset> {
  size_t operator()(const Bitset &x) const { return x.Hash(); }
};
template <>
struct equal_to<Bitset> {
  bool operator()(const Bitset &lhs, const Bitset &rhs) const { return lhs == rhs; }
};
}  // namespace std

// Returns a new Bitset with size equal to `idx_table`'s size. Each entry
// of the new bitset is determined as follows:
// * If the `idx_table` entry is an integer i, then the value is the ith
//   entry of `bitset`.
// * If the entry is nullopt, then the value is False (i.e. zero).
Bitset Remap(const Bitset &bitset, const SizeOptionVector &idx_table);

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Bitset") {
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

  CHECK_EQ(Bitset("1001").ToIndexSetString(), "0,3");
  CHECK_EQ(Bitset("0000").ToIndexSetString(), "");

  auto p = Bitset("000111");
  CHECK_EQ(p.SubsplitChunk(0), Bitset("000"));
  CHECK_EQ(p.SubsplitChunk(1), Bitset("111"));
  CHECK_EQ(p.PCSPChunk(0), Bitset("00"));
  CHECK_EQ(p.PCSPChunk(1), Bitset("01"));
  CHECK_EQ(p.PCSPChunk(2), Bitset("11"));

  CHECK_EQ(Bitset("10011100").RotateSubsplit(), Bitset("11001001"));
  CHECK_EQ(Bitset("010101").SubsplitToIndexSetString(), "1|0,2");

  CHECK_EQ(Bitset("011101").PCSPIsValid(), false);
  CHECK_EQ(Bitset("000111").PCSPIsValid(), false);
  CHECK_EQ(Bitset("100100").PCSPIsValid(), false);
  CHECK_EQ(Bitset("100011001").PCSPIsValid(), true);

  CHECK_EQ(Bitset("100011001").PCSPIsFake(), false);
  CHECK_EQ(Bitset("100011000").PCSPIsFake(), true);

  CHECK_EQ(Bitset("100011001").PCSPParent(), Bitset("100011"));
  CHECK_EQ(Bitset("100011001").PCSPWithoutParent(), Bitset("011001"));
  CHECK_EQ(Bitset("100011001").PCSPChildSubsplit(), Bitset("010001"));
  CHECK_EQ(Bitset("100001110001").PCSPChildSubsplit(), Bitset("01100001"));
  CHECK_EQ(Bitset("100001110001").PCSPChildSubsplitTaxonCounts(), SizePair({1, 2}));
  CHECK_EQ(Bitset("100000111100101").PCSPChildSubsplitTaxonCounts(), SizePair({2, 2}));

  CHECK_EQ(Bitset::Singleton(4, 2), Bitset("0010"));

  // parent clade is 1110, child is 0100, so child subsplit is 1010|0100.
  CHECK_EQ(Bitset::ChildSubsplit(Bitset("00011110"), Bitset("0100")),
           Bitset("10100100"));
  CHECK_EQ(Bitset::ChildSubsplit(Bitset("00011110"), Bitset("1010")),
           Bitset("01001010"));

  CHECK_EQ(Bitset("000110010"), Bitset::PCSPOfPair(Bitset("000110"), Bitset("010100")));
  CHECK_EQ(Bitset("001110010"), Bitset::PCSPOfPair(Bitset("001110"), Bitset("100010")));
  // Missing a left child.
  CHECK_THROWS(Bitset::PCSPOfPair(Bitset("000110"), Bitset("000010")));
  // Missing a right child.
  CHECK_THROWS(Bitset::PCSPOfPair(Bitset("000110"), Bitset("100000")));
  // Not a disjoint union.
  CHECK_THROWS(Bitset::PCSPOfPair(Bitset("000110"), Bitset("100110")));
  // Not a union at all.
  CHECK_THROWS(Bitset::PCSPOfPair(Bitset("000110"), Bitset("100001")));

  CHECK_EQ(true, Bitset("010000").SubsplitIsFake());
  CHECK_EQ(false, Bitset("010010").SubsplitIsFake());
  CHECK_EQ(true, Bitset("").SubsplitIsFake());
  CHECK_EQ(Bitset::FakeSubsplit(Bitset("010")), Bitset("010000"));
  CHECK_EQ(Bitset::FakeChildSubsplit(Bitset("100001")), Bitset("001000"));
  CHECK_THROWS(Bitset::FakeChildSubsplit(Bitset("100011")));
  CHECK_EQ(Bitset::FakePCSP(Bitset("100001")), Bitset("100001000"));
  CHECK_THROWS(Bitset::FakePCSP(Bitset("0000110")));
  CHECK_THROWS(Bitset::FakePCSP(Bitset("100101")));

  // Restrict a bitset.
  CHECK_EQ(Remap(Bitset("10101010101"), {0, 2, 4, 6, 8, 10}), Bitset("111111"));
  // If we apply this remap 3 times we should get back to where we started.
  SizeOptionVector rotate120{6, 7, 8, 0, 1, 2, 3, 4, 5};
  auto to_rotate = Bitset("110010100");
  CHECK_EQ(Remap(Remap(Remap(to_rotate, rotate120), rotate120), rotate120), to_rotate);
  // "Lift" a bitset.
  CHECK_EQ(Remap(Bitset("11"), {0, std::nullopt, 1}), Bitset("101"));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BITSET_HPP_
