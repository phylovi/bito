#include "bitset.hpp"
#include <iostream>
#include <ostream>
#include <sstream>
#include "doctest.h"


Bitset::Bitset(void) : num_set_bits(0) {}


Bitset::Bitset(size_t n, bool def) : value_(n, def), num_set_bits(0) {}


bool Bitset::operator[](size_t i) const { return value_[i]; }


/** Equals comparison */
bool Bitset::operator==(const Bitset& x) const { return x.value_ == value_; }

/** Not-Equals comparison */
bool Bitset::operator!=(const Bitset& x) const {
  return operator==(x) == false;
}


/** Less than comparison */
bool Bitset::operator<(const Bitset& x) const { return x.value_ < value_; }

/** Less than or equals to comparison */
bool Bitset::operator<=(const Bitset& x) const { return x.value_ <= value_; }

/** Greater than comparison */
bool Bitset::operator>(const Bitset& x) const { return x.value_ > value_; }

/** Greater than or equals to comparison */
bool Bitset::operator>=(const Bitset& x) const { return x.value_ < value_; }

/** Bitwise and */
Bitset Bitset::operator&(const Bitset& x) const {
  if (x.value_.size() != value_.size()) {
    throw "Cannot and Bitsets of unequal size";
    abort();
  }
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] && x.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

/** Bitwise or */
Bitset Bitset::operator|(const Bitset& x) const {
  if (x.value_.size() != value_.size()) {
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

/** Bitwise xor */
Bitset Bitset::operator^(const Bitset& x) const {
  if (x.value_.size() != value_.size()) {
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

/** Unary not */
Bitset& Bitset::operator~() {
  for (size_t i = 0; i < value_.size(); i++) {
    value_[i] = !value_[i];
  }

  num_set_bits = value_.size() - num_set_bits;

  return *this;
}

/** Bitwise and assignment */
Bitset& Bitset::operator&=(const Bitset& x) {
  if (x.value_.size() != value_.size()) {
    throw "Cannot and Bitsets of unequal size";
  }

  *this = *this & x;

  return *this;
}

/** Bitwise or assignment */
Bitset& Bitset::operator|=(const Bitset& x) {
  if (x.value_.size() != value_.size()) {
    throw "Cannot or Bitsets of unequal size";
  }

  *this = *this | x;

  return *this;
}


void Bitset::clear(void) {
  // reset the bitset
  value_ = std::vector<bool>(value_.size(), false);
  num_set_bits = 0;
}

bool Bitset::empty(void) const { return value_.empty(); }

void Bitset::flip(size_t i) {
  value_[i] = (value_[i] == false);
  if (value_[i])
    num_set_bits++;
  else
    num_set_bits--;
}

void Bitset::flip() {
  for (size_t i = 0; i < value_.size(); i++) {
    value_[i] = (value_[i] == false);
  }
  // TODO: test
  num_set_bits = value_.size() - num_set_bits;
}

size_t Bitset::getFirstSetBit(void) const {
  size_t index = 0;
  while (index < value_.size() && value_[index] == false) {
    ++index;
  }

  return index;
}

size_t Bitset::getNumberSetBits(void) const {
  // get the internal value_
  return num_set_bits;
}


bool Bitset::isSet(size_t i) const {
  // get the internal value_
  return value_[i];
}

void Bitset::resize(size_t size) { value_.resize(size); }

void Bitset::set(size_t i) {
  if (i >= value_.size()) {
    std::ostringstream ss;
    ss << i;
    throw("Index " + ss.str() +
          " out of bounds in bitset. This will likely cause "
          "unexpected behavior.");
  }

  if (value_[i] == false) {
    ++num_set_bits;
  }

  // set the internal value_
  value_[i] = true;
}


size_t Bitset::size(void) const {
  // get the size from the actual bitset
  return value_.size();
}


void Bitset::unset(size_t i) {
  // TODO(erick): no bounds checking for this one?
  if (value_[i] == true) {
    --num_set_bits;
  }

  value_[i] = false;
}


std::ostream& operator<<(std::ostream& o, const Bitset& x) {
  o << "[";
  for (size_t i = 0; i < x.size(); ++i) {
    o << (x.isSet(i) ? "1" : "0");
  }
  o << "]";

  return o;
}


TEST_CASE("Bitset") {
  auto b = Bitset();
  std::cout << b << std::endl;
}
