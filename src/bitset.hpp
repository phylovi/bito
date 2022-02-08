// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
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

#pragma once

#include <algorithm>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "sugar.hpp"

class Bitset {
 public:
  using BitsetPair = std::pair<Bitset, Bitset>;
  explicit Bitset(std::vector<bool> value);
  // Fills entire Bitset of size `n` with `initial_value`.
  explicit Bitset(size_t n, bool initial_value = false);
  // Builds Bitset from string of "1" and "0"s.
  explicit Bitset(std::string bits_as_str);
  // Builds Bitset of size `n` with only indices in `bits_on` vector set to true.
  explicit Bitset(SizeVector bits_on, size_t n);

  // ** std::bitset Interface Methods
  // These methods are modeled after the std::bitset interface.

  // Get ith bit in Bitset.
  bool operator[](size_t i) const;
  // Get the number of bits in Bitset.
  size_t size() const;
  // Set given index to true.
  void set(size_t i, bool value = true);
  // Sets given index to false.
  void reset(size_t i);
  // Changes each bit to its complement.
  void flip();
  // Comparator:
  // Bitsets are sorted with respect to their binary representation.
  // (e.g. "010" < "101")
  // Note: There are alternative comparators in the Bitset class.
  static int Compare(const Bitset &bitset_a, const Bitset &bitset_b);
  int Compare(const Bitset &other) const;
  // All comparator operator behavior is consistent with std::bitset.
  bool operator==(const Bitset &x) const;
  bool operator!=(const Bitset &x) const;
  bool operator<(const Bitset &x) const;
  bool operator<=(const Bitset &x) const;
  bool operator>(const Bitset &x) const;
  bool operator>=(const Bitset &x) const;
  // All bitwise operator behavior is consistent with std::bitset.
  Bitset operator&(const Bitset &x) const;
  Bitset operator|(const Bitset &x) const;
  Bitset operator^(const Bitset &x) const;
  Bitset operator~() const;
  Bitset operator+(const Bitset &x) const;
  void operator&=(const Bitset &other);
  void operator|=(const Bitset &other);
  // Outputs Bitset string representation to stream.
  friend std::ostream &operator<<(std::ostream &os, const Bitset &bitset);

  // ** Bitset Methods
  // These methods are not from the std::bitset interface.

  // Special Constructors:
  // Make a bitset with only the specified entry turned on.
  static Bitset Singleton(size_t n, size_t which_on);
  // Sets entire bitset to false.
  void Zero();

  // Generates hash value from bitset.
  size_t Hash() const;
  // Outputs bitset as a string of "1" and "0"s.
  std::string ToString() const;
  // Outputs vector of all bit indices set to true.
  SizeVector ToVectorOfSetBits() const;
  // Are all of the bits 1?
  bool All() const;
  // Are any of the bits 1?
  bool Any() const;
  // Are all of the bits 0?
  bool None() const;
  // Get the number of 1s.
  size_t Count() const;
  // Is exactly one of the bits 1?
  bool IsSingleton() const;
  // Is this disjoint with the given bitset?
  bool IsDisjoint(const Bitset &other) const;
  // Take the minimum of the bitset and its complement.
  void Minorize();
  // Copy all of the bits from another bitset into this bitset, starting at
  // begin, and optionally flipping the bits as they get copied.
  void CopyFrom(const Bitset &other, size_t begin, bool flip);
  // If the bitset only has one bit on, then we return the location of that bit.
  // Otherwise, return nullopt.
  std::optional<uint32_t> SingletonOption() const;
  // Output as string of comma-separated indices.
  std::string ToVectorOfSetBitsAsString() const;
  // Returns a new Bitset with size equal to `idx_table`'s size. Each entry
  // of the new bitset is determined as follows:
  // - If the `idx_table` entry is an integer i, then the value is the ith
  //   entry of `bitset`.
  // - If the entry is nullopt, then the value is False (i.e. zero).
  static Bitset Remap(const Bitset &bitset, const SizeOptionVector &idx_table);

  // ** Clade / MultiClade Methods
  // These methods require bitsets to represent "clades". A clade is an expression of a
  // subset of a taxon set.  The size of the total taxon set is equal to the size of the
  // clade's bitset, with each bit index representing a inclusion/exclusion of a
  // specific member of that taxon set.
  //
  // There are bitset "types" composed of multiple clades (here, called a
  // MultiClade). Subsplits are composed of two clades and Edges are composed of three
  // clades.

  // Comparator: Clades are sorted with respect to the lexigraphical representation of
  // their taxon subset. (e.g. If two clades are "010" and "101", then their taxon
  // subsets are {b} and {a,c}. "b" > "a", therefore "010" > "101".) Note: Sorting by
  // taxon representation gives the precise opposite ordering to sorting by binary
  // representation (except in the case of zero/empty set!).
  static int CladeCompare(const Bitset &bitset_a, const Bitset &bitset_b);
  int CladeCompare(const Bitset &other) const;
  // For a specified number of clades, return the clade/taxon size.
  size_t MultiCladeGetCladeSize(const size_t clade_count) const;
  // For a specified number of clades, return a bitset of the ith clade.
  Bitset MultiCladeGetClade(const size_t which_clade, const size_t clade_count) const;
  // For a specified number of clades, return a string of "1" and "0" for each clade,
  // with clades separated by a "|".
  std::string MultiCladeToString(const size_t clade_count) const;

  // ** Subsplit Methods
  // These methods require bitset to represent "subsplits".  A subsplit is of even
  // length, consisting of two equal-sized, disjoint "clades", representing the two
  // sides of the subsplit. Clades are normally stored in a sorted order wrt to their
  // lexicographic taxon ordering: the smaller "left" (or "sorted") clade stored in the
  // 0-position, and the larger "right" (or "rotated") clade in the 1-position.

  static inline size_t SubsplitCladeCount = 2;
  enum class SubsplitClade : bool {
    Left = 0,  
    Right = 1, 
  };

  // Constructors:
  // Each argument represents one of the clade of the subsplit.
  // Each clade follows the corresponding Bitset constructor.
  //
  // Build a Subsplit bitset out of a compatible pair of clades.
  static Bitset Subsplit(const Bitset &clade_0, const Bitset &clade_1);
  // Builds Clades from strings of "1" and "0" characters.
  static Bitset Subsplit(const std::string clade_0, const std::string clade_1);
  // Builds Clades of size `n` with only indices in the clade vectors set to true.
  static Bitset Subsplit(const SizeVector clade_0, const SizeVector clade_1,
                         const size_t n);
  // Given two arbitrarily ordered clades of same length, return a subsplit with the
  // clades in sorted order by taxon representation.
  static Bitset SubsplitFromUnorderedClades(const Bitset &clade_0,
                                            const Bitset &clade_1);
  // Special Constructors:
  // A "leaf subsplit" is a subsplit where one of the clades is equal
  // to the empty set (zero). In practice, these are allowed in only two cases: When
  // subsplit is a leaf or root.  A leaf should only have a single member in its
  // non-empty clade, and a root should have the full taxon set in its non-empty clade.
  //
  // Make a "leaf" subsplit (pairs given nonempty_clade with an empty_clade).
  static Bitset LeafSubsplitOfNonemptyClade(const Bitset &nonempty_clade);
  // Make a "leaf" child subsplit of a given parent subsplit. The leftside clade of
  // the parent subsplit must be non-empty and the rightside clade must be a singleton.
  static Bitset LeafSubsplitOfParentSubsplit(const Bitset &parent_subsplit);
  // Make the UCA (universal common ancestor) subsplit of the DAG root node with the given taxon count.
  // Since subsplit bitsets are always big-small, the DAG root node subsplit
  // consists of all 1s then 0s (e.g. 5 would return '11111|00000').
  static Bitset UCASubsplitOfTaxonCount(const size_t taxon_count);
  // Get the full rootsplit bitset out of a rootsplit half.
  // Note: the first half of the rootsplit bitset is always larger than the second.
  static Bitset RootsplitOfClade(const Bitset &clade);
  // Comparator:
  // Subsplits are sorting on the following:
  // (1) The number of taxa in each of their subsplits.
  // (2) The std::bitset ordering of each of their respective unions.
  // (3) The std::bitset ordering of each or their sorted clades.
  static int SubsplitCompare(const Bitset &subsplit_a, const Bitset &subsplit_b);
  int SubsplitCompare(const Bitset &other) const;
  // Flip the order of the two clades of a subsplit.
  Bitset SubsplitRotate() const;
  // Sorts clades of subsplit so that they are ordered by their taxon representation.
  Bitset SubsplitSortClades() const;
  // Gets the size of each of each clade. This is the same as the size of the whole
  // taxon set.
  size_t SubsplitGetCladeSize() const;
  // Get clade according to its taxon ordering.
  Bitset SubsplitGetClade(const size_t which_clade) const;
  Bitset SubsplitGetClade(const SubsplitClade which_clade) const;
  // Output subsplit as string of "1" and "0" characters, with each clade separated by a
  // "|".
  std::string SubsplitToString() const;
  // Output subsplit to string as a comma-separated list of true bits positions, with
  // each clade separated by a "|".
  std::string SubsplitToVectorOfSetBitsAsString() const;
  // Is this the subsplit of a leaf node?
  bool SubsplitIsLeaf() const;
  // Is this the subsplit of root node?
  bool SubsplitIsRoot() const;
  // Is this the subsplit of a rootsplit?
  bool SubsplitIsRootsplit() const;
  // Is this the left/rotated clade of the given subsplit?
  bool SubsplitIsLeftChildOf(const Bitset &parent) const;
  // Is this the right/sorted clade of the given subsplit?
  bool SubsplitIsRightChildOf(const Bitset &parent) const;
  // Get the union of the two clades.
  Bitset SubsplitCladeUnion() const;
  // Get whether the given child is the sorted (false, 0-position) or rotated (true,
  // 1-position) child to the given parent.
  static bool SubsplitIsWhichChildOf(const Bitset &parent, const Bitset &child);
  // Check whether subsplits form a child/parent pair.
  static bool SubsplitIsParentChildPair(const Bitset &parent, const Bitset &child);
  // Check whether subsplits are adjacent/related (whether either is the parent of the
  // other).
  static bool SubsplitIsAdjacent(const Bitset &subsplit_a, const Bitset &subsplit_b);
  // Check whether bitset represents valid Subsplit (contains two equal-sized, disjoint
  // clades).
  bool SubsplitIsValid() const;

  // ** Edge methods
  // These functions require the bitset to be a "edge bitset"  (also called PCSP, or
  // parent-child subsplit pair). Edges are composed of three equal-sized "clades". The
  // first clade represents the sister clade, the second the focal clade, and the third
  // clade describes the "sorted" clade of the child subsplit. We define "sorted" as the
  // clade of a child subsplit that has a bitset with the smaller binary representation.
  // The clades are well defined relative to the cut parent: the other part of
  // the subsplit, "rotated clade" is the cut parent set, minus the "sorted clade".
  //
  // For example, `100|011|001` is composed of the clades `100`, `011` and `001`.
  // If the taxa are x0, x1, and x2 then this means the parent subsplit is (A,
  // BC) with bitset `100|011`, and the child subsplit is (B, C) with bitset
  // `010|001.` Child_0 is the clade `001` and child_1 is the clade `010.`
  //
  // For rootsplit Edges where the parent subsplit is the DAG root node, the Edge
  // is the sister clade (all 0s), the focal clade (all 1s), and "clade 0". For
  // example, `000111010` is the Edge from the DAG root node to the rootsplit (AC, B).
  // See the unit tests at the bottom for more examples.

  static inline size_t PCSPCladeCount = 3;
  enum class PCSPClade : size_t {
    Sister = 0,
    Focal = 1,
    LeftChild = 2
  };

  // Constructors:
  // Build a Edge bitset from a compatible parent-child pair of
  // Subsplit bitsets.
  static Bitset Edge(const Bitset &parent_subsplit, const Bitset &child_subsplit);
  // Build a Edge bitset from explicit sister, focal, and sorted-child clades.
  static Bitset Edge(const Bitset &sister_clade, const Bitset &focal_clade,
                     const Bitset &sorted_child_clade);
  // Builds sister, focal, and sorted-child clades from strings of "1" and "0"
  // characters.
  static Bitset Edge(const std::string sister_clade, const std::string focal_clade,
                     const std::string sorted_child_clade);
  // Special Constructors:
  // Make a "leaf" Edge of a given parent subsplit; assert that the left-hand clade of
  // the parent subsplit is non-empty and that the right-hand clade is a singleton.
  // This leaf subsplit has parent subsplit on the left and all zeroes on the right.
  static Bitset EdgeToLeaf(const Bitset &parent_subsplit);
  // Given a rootsplit, get the Edge connecting the DAG root node to that rootsplit
  // (e.g. '1100|0011' would return '0000|1111|0011').
  static Bitset EdgeToRootsplit(const Bitset &rootsplit);
  // Output Edge as string of "1" and "0" characters, with each clade separated by a
  // "|".
  std::string EdgeToString() const;
  // Checks whether bitset represents a valid set of taxon clades for Edge.
  bool EdgeIsValid() const;
  // Checks whether the Edge clade-side edge goes to a leaf subsplit.
  bool EdgeIsLeaf() const;
  // Sorts Edge so that parent and child are rotated properly so that second clade is
  // the focal clade of the parent and third clade is the sorted side of the child.
  Bitset EdgeSortClades() const;
  // Do the sister and focal clades union to the whole taxon set?
  // Method excludes rootsplit Edges where sister and focal clades
  // also union to the whole taxon set.
  bool EdgeIsParentRootsplit() const;
  // Gets the size of each of each clade. This is the same as the size of the whole
  // taxon set.
  size_t EdgeGetCladeSize() const;
  // Get the ith clade of the Edge.
  Bitset EdgeGetClade(const size_t which_clade) const;
  Bitset EdgeGetClade(const PCSPClade which_clade) const;
  // Get the parent subsplit of the Edge.
  Bitset EdgeGetParentSubsplit() const;
  // Get the child subsplit of the Edge.
  Bitset EdgeGetChildSubsplit() const;
  // Get the number of taxa in each side of the child subsplit.
  SizePair EdgeGetChildSubsplitTaxonCounts() const;

 protected:
  // Vector of bits.
  std::vector<bool> value_;
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

TEST_CASE("Bitset: Clades, Subsplits, Edges") {
  auto p = Bitset("000111");
  // Subsplit: 000|111
  CHECK_EQ(p.SubsplitGetClade(Bitset::SubsplitClade::Left), Bitset("000"));
  CHECK_EQ(p.SubsplitGetClade(Bitset::SubsplitClade::Right), Bitset("111"));
  // Edge: 00|01|11
  CHECK_EQ(p.EdgeGetClade(Bitset::PCSPClade::Sister), Bitset("00"));
  CHECK_EQ(p.EdgeGetClade(Bitset::PCSPClade::Focal), Bitset("01"));
  CHECK_EQ(p.EdgeGetClade(Bitset::PCSPClade::LeftChild), Bitset("11"));

  CHECK_EQ(Bitset("11001010").SubsplitCladeUnion(), Bitset("1110"));

  CHECK_EQ(Bitset("10011100").SubsplitRotate(), Bitset("11001001"));
  CHECK_EQ(Bitset("010101").SubsplitToVectorOfSetBitsAsString(), "1|0,2");

  CHECK_EQ(Bitset("101010").SubsplitIsLeftChildOf(Bitset("111000")), true);
  // CHECK_EQ(Bitset::SubsplitIsWhichChildOf(Bitset("111000"), Bitset("101010")),true);
  CHECK_EQ(Bitset("00100001").SubsplitIsRightChildOf(Bitset("11000011")), true);
  // CHECK_EQ(Bitset::SubsplitIsWhichChildOf(Bitset("000111"), Bitset("101010")),
  // false);
  CHECK_EQ(Bitset("010001").SubsplitIsLeftChildOf(Bitset("110001")), false);
  CHECK_EQ(Bitset("010001").SubsplitIsRightChildOf(Bitset("01000011")), false);
  // Should throw because Bitsets can't be divided into equal-sized clades.
  CHECK_THROWS(Bitset("11010").SubsplitIsLeftChildOf(Bitset("10101")));
  CHECK_THROWS(Bitset("11010").SubsplitIsRightChildOf(Bitset("10101")));

  CHECK_EQ(Bitset("101010").SubsplitIsRootsplit(), true);
  CHECK_EQ(Bitset("111000").SubsplitIsRootsplit(), false);
  CHECK_EQ(Bitset("11000001").SubsplitIsRootsplit(), false);

  CHECK_EQ(Bitset("011101").EdgeIsValid(), false);
  CHECK_EQ(Bitset("000111").EdgeIsValid(), false);
  CHECK_EQ(Bitset("100100").EdgeIsValid(), false);
  CHECK_EQ(Bitset("100011001").EdgeIsValid(), true);

  CHECK_EQ(Bitset("100011001").EdgeIsLeaf(), false);
  CHECK_EQ(Bitset("100011000").EdgeIsLeaf(), true);

  CHECK_EQ(Bitset("000111010").EdgeIsParentRootsplit(), false);
  CHECK_EQ(Bitset("000111000100").EdgeIsParentRootsplit(), false);
  CHECK_EQ(Bitset("101010000").EdgeIsParentRootsplit(), true);

  CHECK_EQ(Bitset("100011001").EdgeGetParentSubsplit(), Bitset("100011"));
  CHECK_EQ(Bitset("011100001").EdgeGetParentSubsplit(), Bitset("100011"));
  CHECK_EQ(Bitset("100011001").EdgeGetChildSubsplit(), Bitset("010001"));
  CHECK_EQ(Bitset("100001110001").EdgeGetChildSubsplit(), Bitset("01100001"));
  CHECK_EQ(Bitset("100001110001").EdgeGetChildSubsplitTaxonCounts(), SizePair({1, 2}));
  CHECK_EQ(Bitset("100000111100101").EdgeGetChildSubsplitTaxonCounts(),
           SizePair({2, 2}));

  CHECK_EQ(Bitset::Singleton(4, 2), Bitset("0010"));

  CHECK_EQ(Bitset("100010"), Bitset::Subsplit(Bitset("100"), Bitset("010")));
  CHECK_EQ(Bitset("110001"), Bitset::Subsplit(Bitset("001"), Bitset("110")));
  // Invalid clade pair.
  CHECK_THROWS(Bitset::Subsplit(Bitset("1100"), Bitset("001")));
  CHECK_THROWS(Bitset::Subsplit(Bitset("111"), Bitset("001")));

  CHECK_EQ(Bitset("000110010"), Bitset::Edge(Bitset("110000"), Bitset("100010")));
  CHECK_EQ(Bitset("110001000"), Bitset::Edge(Bitset("110001"), Bitset("001000")));
  // Invalid parent-child pair.
  CHECK_THROWS(Bitset::Edge(Bitset("110001"), Bitset("010001")));
  CHECK_THROWS(Bitset::Edge(Bitset("11000101"), Bitset("010001")));
  CHECK_THROWS(Bitset::Edge(Bitset("110001"), Bitset("110100")));

  CHECK_EQ(Bitset::RootsplitOfClade(Bitset("0011")), Bitset("11000011"));
  CHECK_EQ(Bitset::EdgeToRootsplit(Bitset("11000011")), Bitset("000011110011"));

  CHECK_EQ(Bitset("010000").SubsplitIsLeaf(), true);
  CHECK_EQ(Bitset("010010").SubsplitIsLeaf(), false);
  CHECK_EQ(Bitset("111000").SubsplitIsLeaf(), false);
  CHECK_EQ(Bitset::LeafSubsplitOfNonemptyClade(Bitset("010")), Bitset("010000"));
  CHECK_EQ(Bitset::LeafSubsplitOfParentSubsplit(Bitset("100001")), Bitset("001000"));
  CHECK_THROWS(Bitset::LeafSubsplitOfParentSubsplit(Bitset("100011")));
  CHECK_EQ(Bitset::EdgeToLeaf(Bitset("100001")), Bitset("100001000"));
  CHECK_THROWS(Bitset::EdgeToLeaf(Bitset("0000110")));
  CHECK_THROWS(Bitset::EdgeToLeaf(Bitset("100101")));

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
