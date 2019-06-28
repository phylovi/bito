#include "bitset.hpp"
#include <cassert>


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

size_t Bitset::size(void) const { return value_.size(); }

void Bitset::set(size_t i, bool value) {
  assert(i < value_.size());
  value_[i] = value;
}

void Bitset::reset(size_t i) {
  assert(i < value_.size());
  value_[i] = false;
}

void Bitset::flip() { value_.flip(); }

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
  assert(value_.size() == x.size());
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] && x.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

Bitset Bitset::operator|(const Bitset& x) const {
  assert(value_.size() == x.size());
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] || x.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

Bitset Bitset::operator^(const Bitset& x) const {
  assert(value_.size() == x.size());
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] != x.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

Bitset Bitset::operator~() const {
  Bitset r(value_);
  r.value_.flip();
  return r;
}

void Bitset::operator&=(const Bitset& other) {
  assert(value_.size() == other.size());
  for (size_t i = 0; i < value_.size(); i++) {
    value_[i] = value_[i] && other[i];
    }
}

void Bitset::operator|=(const Bitset& other) {
  assert(value_.size() == other.size());
  for (size_t i = 0; i < value_.size(); i++) {
    value_[i] = value_[i] || other[i];
  }
}

// These methods aren't in the bitset interface.

size_t Bitset::Hash(void) const {
  return std::hash<std::vector<bool>>{}(value_);
}

std::string Bitset::ToString() const {
  std::string str;
  for (size_t i = 0; i < value_.size(); ++i) {
    str += (value_[i] ? '1' : '0');
  }
  return str;
}

void Bitset::Minorize() {
  assert(value_.size() > 0);
  if (value_[0]) {
    value_.flip();
  }
}


