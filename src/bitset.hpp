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
  // Builds Bitset from boolean vector.
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

  // Get underlying data vector.
  std::vector<bool> GetData();

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
  // MultiClade). Subsplits are composed of two clades and PCSPs are composed of three
  // clades.

  // Comparator: Clades are ordered with respect to the lexigraphical representation of
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
  // These methods require bitset to represent "subsplits".  Subsplits represent nodes
  // within the SubsplitDAG.  A subsplit are composed of two equal-sized, disjoint
  // "clades", representing a fork in the DAG and the taxon sets descending from the
  // left and right sides of the fork. Clades are normally stored in a sorted order wrt
  // to their lexicographic taxon ordering: the smaller "left" clade stored in the
  // 0-position, and the larger "right" clade in the 1-position.

  enum class SubsplitClade : size_t { Left, Right, Unspecified };
  static const inline size_t SubsplitCladeCount = 2;
  class SubsplitCladeEnum
      : public EnumWrapper<SubsplitClade, size_t, SubsplitCladeCount,
                           SubsplitClade::Left, SubsplitClade::Right> {
   public:
    static inline const std::string Prefix = "SubsplitClade";
    static inline const Array<std::string> Labels = {{"Left", "Right"}};

    static std::string ToString(const SubsplitClade e) {
      std::stringstream ss;
      ss << Prefix << "::" << Labels[e];
      return ss.str();
    }
    friend std::ostream &operator<<(std::ostream &os, const SubsplitClade e) {
      os << ToString(e);
      return os;
    }
  };
  // Does not iterate over "unspecified" clade.
  using SubsplitCladeIterator =
      EnumIterator<SubsplitClade, SubsplitClade::Left, SubsplitClade::Right>;

  static SubsplitClade Opposite(const SubsplitClade clade) {
    switch (clade) {
      case SubsplitClade::Left:
        return SubsplitClade::Right;
      case SubsplitClade::Right:
        return SubsplitClade::Left;
      default:
        Failwith("Cannot get Opposite of Unspecified Clade");
    }
  }

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
  // Make a "leaf subsplit" (pairs given nonempty_clade with an empty_clade).
  static Bitset LeafSubsplitOfNonemptyClade(const Bitset &nonempty_clade);
  // Make a "leaf subsplit" of a given parent subsplit. The left clade of
  // the parent subsplit must be non-empty and the right clade must be a singleton.
  static Bitset LeafSubsplitOfParentSubsplit(const Bitset &parent_subsplit);
  // Make the UCA (universal common ancestor) subsplit of the DAG root node with the
  // given taxon count. Since subsplit bitsets are always big-small, the DAG root node
  // subsplit consists of all 1s then 0s (e.g. 5 would return '11111|00000').
  static Bitset UCASubsplitOfTaxonCount(const size_t taxon_count);

  // Get the full rootsplit bitset out of a rootsplit half.
  // Note: the first half of the rootsplit bitset is always larger than the second.
  static Bitset RootsplitSubsplitOfClade(const Bitset &clade);
  // Comparator:
  // Subsplits are sorting on the following:
  // (1) The number of taxa in each of their subsplits.
  // (2) The std::bitset ordering of each of their respective unions.
  // (3) The std::bitset ordering of each of their left clades.
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
  // #350 why do we want these two overloads?
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
  // Is this the UCA subsplit?
  bool SubsplitIsUCA() const;
  // Is this the subsplit of a rootsplit?
  bool SubsplitIsRootsplit() const;
  // Is this the left clade of the given subsplit?
  bool SubsplitIsLeftChildOf(const Bitset &parent) const;
  // Is this the right clade of the given subsplit?
  bool SubsplitIsRightChildOf(const Bitset &parent) const;
  // Get the union of the two clades.
  Bitset SubsplitCladeUnion() const;
  // Get whether the given child is the left or right child to the given parent.
  static SubsplitClade SubsplitIsChildOfWhichParentClade(const Bitset &parent,
                                                         const Bitset &child);
  // Check whether subsplits form a child/parent pair.
  static bool SubsplitIsParentChildPair(const Bitset &parent, const Bitset &child);
  // Check whether subsplits are adjacent/related (whether either is the parent of the
  // other).
  static bool SubsplitIsAdjacent(const Bitset &subsplit_a, const Bitset &subsplit_b);
  // Check whether subsplits are potentially related (ancestor/descendant) pair (union
  // of descendant is a subset of one of clades of ancestor).
  static bool SubsplitIsAncestorDescendantPair(const Bitset &ancestor,
                                               const Bitset &descendant,
                                               const SubsplitClade clade_type);
  // Check whether bitset represents valid Subsplit (contains two equal-sized,
  // disjoint clades).
  bool SubsplitIsValid() const;

  // ** PCSP methods
  // These functions require the bitset to be a "PCSP bitset" (parent-child subsplit
  // pair). PCSP represent edges between nodes within the SubsplitDAG. They are composed
  // of three equal-sized "clades": (0) sister clade of parent, (1) focal clade of
  // parent, (2) the right clade of the child. We define the "right" clade of a child
  // subsplit that has a bitset with the larger lexicographic representation. The
  // remaining clade are well-defined relative to the focal parent subsplit.
  //
  // For example, `100|011|001` is composed of the clades `100`, `011` and `001`.
  // If the taxa are x0, x1, and x2 then this means the parent subsplit is (A,
  // BC) with bitset `100|011`, and the child subsplit is (B, C) with bitset
  // `010|001.` Child_0 is the clade `001` and child_1 is the clade `010.`
  //
  // For rootsplit PCSPs where the parent subsplit is the DAG root node, the PCSP
  // is the sister clade (all 0s), the focal clade (all 1s), and "clade 0". For
  // example, `000111010` is the PCSP from the DAG root node to the rootsplit (AC, B).
  // See the unit tests at the bottom for more examples.

  static inline size_t PCSPCladeCount = 3;
  enum class PCSPClade : size_t { Sister, Focal, RightChild };
  using PCSPCladeIterator =
      EnumIterator<PCSPClade, PCSPClade::Sister, PCSPClade::RightChild>;

  // Constructors:
  // Build a PCSP bitset from a compatible parent-child pair of
  // Subsplit bitsets.
  static Bitset PCSP(const Bitset &parent_subsplit, const Bitset &child_subsplit);
  // Build a PCSP bitset from explicit sister, focal, and right child clades.
  static Bitset PCSP(const Bitset &sister_clade, const Bitset &focal_clade,
                     const Bitset &right_child_clade);
  // Builds sister, focal, and right child clades from strings of "1" and "0"
  // characters.
  static Bitset PCSP(const std::string sister_clade, const std::string focal_clade,
                     const std::string right_child_clade);
  // Special Constructors:
  // Make a PCSP from parent subsplit to child leaf subsplit. Asserts that the left-hand
  // clade of the parent subsplit is non-empty and that the right-hand clade is a
  // singleton. This leaf subsplit has parent subsplit on the left and all zeroes on the
  // right.
  static Bitset PCSPFromRightParentCladeToLeaf(const Bitset &parent_subsplit);
  // Given a rootsplit, get the PCSP connecting the DAG root node to that rootsplit
  // (e.g. '1100|0011' would return '0000|1111|0011').
  static Bitset PCSPFromUCAToRootsplit(const Bitset &rootsplit);
  // Output PCSP as string of "1" and "0" characters, with each clade separated by a
  // "|".
  std::string PCSPToString() const;
  // Checks whether bitset represents a valid set of taxon clades for PCSP.
  bool PCSPIsValid() const;
  // Checks whether the PCSP child is a leaf subsplit.
  bool PCSPChildIsLeaf() const;
  // Sorts PCSP so that parent and child are arranged properly so that second clade is
  // the focal clade of the parent and third clade is the right side of the child.
  Bitset PCSPSortClades() const;
  // Do the sister and focal clades union to the whole taxon set?
  // Method excludes rootsplit PCSPs where sister and focal clades
  // also union to the whole taxon set.
  bool PCSPIsParentRootsplit() const;
  // Gets the size of each of each clade. This is the same as the size of the whole
  // taxon set.
  size_t PCSPGetCladeSize() const;
  // Get the ith clade of the PCSP.
  Bitset PCSPGetClade(const size_t which_clade) const;
  Bitset PCSPGetClade(const PCSPClade which_clade) const;
  // Get the parent subsplit of the PCSP.
  Bitset PCSPGetParentSubsplit() const;
  // Get the child subsplit of the PCSP.
  Bitset PCSPGetChildSubsplit() const;
  // Get the number of taxa in each side of the child subsplit.
  SizePair PCSPGetChildSubsplitTaxonCounts() const;
  // Checks if PCSP are adjacent. In other words, PCSPs share a common node.
  static bool PCSPIsParentChildPair(const Bitset &parent_pcsp,
                                    const Bitset &child_pcsp);

 protected:
  // Vector of bits.
  std::vector<bool> value_;
};

using Subsplit = Bitset;
using PCSP = Bitset;

using SubsplitClade = Bitset::SubsplitClade;
using SubsplitCladeEnum = Bitset::SubsplitCladeEnum;
using PCSPClade = Bitset::PCSPClade;

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
