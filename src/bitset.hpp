// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BITSET_HPP_
#define SRC_BITSET_HPP_

#include <algorithm>
#include <functional>
#include <string>
#include <vector>
#include "optional.hpp"

// A rewrite of the RbBitSet class from RevBayes by Sebastian Hoehna.
// In general, I'm trying to follow the interface of std::bitset.

class Bitset {
 public:
  explicit Bitset(std::vector<bool> value);
  explicit Bitset(size_t n, bool initial_value = false);
  explicit Bitset(std::string);

  Bitset copy();
  // // Eliminate copy constructors.
  // Bitset(const Bitset &) = delete;
  // Bitset &operator=(const Bitset &) = delete;
  // // Keep move constructor. This is necessary because if we fiddle with copy
  // // constructors it drops the default copy constructor too.
  // Bitset(Bitset &&) = default;

  bool operator[](size_t i) const;
  size_t size(void) const;

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

  void operator&=(const Bitset &other);
  void operator|=(const Bitset &other);

  // These methods aren't in the bitset interface, so they get our usual
  // convention.
  size_t Hash(void) const;
  std::string ToString() const;
  bool Any() const;
  void Minorize();
  void CopyFrom(const Bitset &other, size_t begin, bool flip);
  // If the bitset only has one bit on, then we return the location of that bit.
  // Otherwise, return nullopt.
  std::experimental::optional<uint32_t> SingletonOption() const;

  Bitset SisterExchange() const;
  Bitset SplitChunk(size_t i) const;
  std::string ToStringChunked(size_t chunk_count) const;
  std::string SubsplitToString() const;
  std::string PCSSToString() const;
  bool PCSSIsValid() const;
  size_t PCSSChunkSize() const;
  Bitset PCSSChunk(size_t i) const;
  Bitset PCSSParent() const;
  Bitset PCSSChild() const;

  // Static methods
  // Make the full subsplit out of the parent subsplit, whose second half is
  // split by child_half.
  static Bitset ChildSubsplit(const Bitset &parent_subsplit,
                              const Bitset &child_half);

 private:
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
  bool operator()(const Bitset &lhs, const Bitset &rhs) const {
    return lhs == rhs;
  }
};
}  // namespace std

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
  CHECK_EQ(build_up, Bitset("0101"));

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
  CHECK_EQ(std::min(Bitset("1100"), Bitset("1010")), Bitset("1010"));

  a &= Bitset("0110");
  CHECK_EQ(a, Bitset("0100"));

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

  auto p = Bitset("000111");
  CHECK_EQ(p.PCSSChunk(0), Bitset("00"));
  CHECK_EQ(p.PCSSChunk(1), Bitset("01"));
  CHECK_EQ(p.PCSSChunk(2), Bitset("11"));

  CHECK_EQ(Bitset("10011100").SisterExchange(), Bitset("11001001"));

  CHECK_EQ(Bitset("011101").PCSSIsValid(), false);
  CHECK_EQ(Bitset("000111").PCSSIsValid(), false);
  CHECK_EQ(Bitset("100100").PCSSIsValid(), false);
  CHECK_EQ(Bitset("100011001").PCSSIsValid(), true);

  CHECK_EQ(Bitset("100011001").PCSSParent(), Bitset("100011"));
  CHECK_EQ(Bitset("100011001").PCSSChild(), Bitset("010001"));
  CHECK_EQ(Bitset("100001110001").PCSSChild(), Bitset("01100001"));

  // parent clade is 1110, child is 0100, so child subsplit is 1010|0100.
  CHECK_EQ(Bitset::ChildSubsplit(Bitset("00011110"), Bitset("0100")),
           Bitset("10100100"));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BITSET_HPP_
