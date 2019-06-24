#include "bitset.hpp"
#include <cassert>
#include "doctest.h"


Bitset::Bitset(void) {}

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

void Bitset::set(size_t i) {
  assert(i < value_.size());
  value_[i] = true;
}

void Bitset::reset(size_t i) {
  assert(i < value_.size());
  value_[i] = false;
}

size_t Bitset::size(void) const { return value_.size(); }

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


TEST_CASE("Bitset") {
  auto a = Bitset("1100");

  CHECK(a[2] == false);
  CHECK(a[1] == true);

  auto build_up = Bitset(4);
  build_up.set(1);
  build_up.set(3);
  CHECK(build_up == Bitset("0101"));

  auto strip_down = Bitset(4, true);
  strip_down.reset(0);
  strip_down.reset(2);
  CHECK(build_up == Bitset("0101"));

  CHECK(a.size() == 4);

  CHECK(Bitset("1100") == Bitset("1100"));
  CHECK(Bitset("1100") != Bitset("0100"));
  CHECK(Bitset("0100") < Bitset("0110"));
  CHECK(Bitset("0100") < Bitset("0110"));
  CHECK(Bitset("0010") < Bitset("0100"));

  CHECK((Bitset("1100") & Bitset("1010")) == Bitset("1000"));
  CHECK((Bitset("1100") | Bitset("1010")) == Bitset("1110"));
  CHECK((Bitset("1100") ^ Bitset("1010")) == Bitset("0110"));
  CHECK(~Bitset("1010") == Bitset("0101"));
}
