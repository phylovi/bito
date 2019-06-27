#include "bitset.hpp"
#include <cassert>
#include "doctest.h"

// A rewrite of the RbBitSet class from RevBayes by Sebastian Hoehna.
// In general, I'm trying to follow the interface of std::bitset.


Bitset::Bitset(std::vector<bool> value) : value_(value) {}

Bitset::Bitset(size_t n, bool initial_value) : value_(n, initial_value) {}

Bitset::Bitset(std::string str) : Bitset(str.length()) {
  for (size_t i = 0; i < value_.size(); i++) {
    if (str[i] == '0') {
      value_[i] = false;
    } else if (str[i] == '1') {
      value_[i] = true;
    } else {
      throw("String constructor for Bitset must use only 0s or 1s; found '" +
            std::string(1, str[i]) + "'.");
    }
  }
}

bool Bitset::operator[](size_t i) const { return value_[i]; }

void Bitset::set(size_t i, bool value) {
  assert(i < value_.size());
  value_[i] = value;
}

void Bitset::reset(size_t i) {
  assert(i < value_.size());
  value_[i] = false;
}

size_t Bitset::size(void) const { return value_.size(); }
size_t Bitset::hash(void) const {
  return std::hash<std::vector<bool>>{}(value_);
}

bool Bitset::operator==(const Bitset& other) const {
  return value_ == other.value_;
}
bool Bitset::operator!=(const Bitset& other) const {
  return value_ != other.value_;
}
bool Bitset::operator<(const Bitset& other) const {
  return value_ < other.value_;
}
bool Bitset::operator<=(const Bitset& other) const {
  return value_ <= other.value_;
}
bool Bitset::operator>(const Bitset& other) const {
  return value_ > other.value_;
}
bool Bitset::operator>=(const Bitset& other) const {
  return value_ >= other.value_;
}

Bitset Bitset::operator&(const Bitset& x) const {
  if (value_.size() != x.size()) {
    throw "Cannot and Bitsets of unequal size";
  }
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] && x.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

// Let x be the and of x and and_with_this.
void Bitset::AndWith(Bitset& x, const Bitset& and_with_this) {
  if (x.size() != and_with_this.size()) {
    throw "Cannot and Bitsets of unequal size";
  }
  for (size_t i = 0; i < x.size(); i++) {
    x.set(i, x[i] && and_with_this[i]);
    }
}

// Bitwise or
Bitset Bitset::operator|(const Bitset& x) const {
  if (value_.size() != x.size()) {
    throw "Cannot or Bitsets of unequal sizes";
  }
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] || x.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

// Bitwise xor
Bitset Bitset::operator^(const Bitset& x) const {
  if (value_.size() != x.size()) {
    throw "Cannot xor Bitsets of unequal size";
  }
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] != x.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

// Unary not
Bitset Bitset::operator~() const {
  Bitset r(value_);
  r.value_.flip();
  return r;
}

std::string Bitset::ToString() {
  std::string str = "[";
  for (size_t i = 0; i < value_.size(); ++i) {
    str += (value_[i] ? '1' : '0');
  }
  str += ']';
  return str;
}

// This is how we inject a hash routine into the std namespace so that we can
// use it as a key for an unordered_map.
// https://en.cppreference.com/w/cpp/container/unordered_map
namespace std {
template <>
struct hash<Bitset> {
  size_t operator()(Bitset const& x) const noexcept { return x.hash(); }
};
}

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
}
