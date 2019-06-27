#ifndef SRC_BITSET_HPP_
#define SRC_BITSET_HPP_

#include <string>
#include <vector>

// A rewrite of the RbBitSet class from RevBayes by Sebastian Hoehna.
// In general, I'm trying to follow the interface of std::bitset.

class Bitset {
 public:
  explicit Bitset(std::vector<bool> value);
  explicit Bitset(size_t n, bool initial_value = false);
  explicit Bitset(std::string);

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
  std::string ToString();
  void Minorize();

 private:
  std::vector<bool> value_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Bitset") {
  auto a = Bitset("1100");

  CHECK_EQ(a[2], false);
  CHECK_EQ(a[1], true);

  auto build_up = Bitset(4);
  build_up.set(1);
  build_up.set(3);
  CHECK_EQ(build_up, Bitset("0101"));

  auto strip_down = Bitset(4, true);
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

  a &= Bitset("0110");
  CHECK_EQ(a, Bitset("0100"));

  a.flip();
  CHECK_EQ(a, Bitset("1011"));
  a.Minorize();
  CHECK_EQ(a, Bitset("0100"));
  a.Minorize();
  CHECK_EQ(a, Bitset("0100"));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BITSET_HPP_
