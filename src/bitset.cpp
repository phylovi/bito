// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

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

// Put here in anticipation of turning off copy constructor.
Bitset Bitset::copy() { return Bitset(this->value_); }

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

Bitset Bitset::operator+(const Bitset& other) const {
  Bitset sum(value_.size() + other.size());
  sum.CopyFrom(*this, 0, false);
  sum.CopyFrom(other, value_.size(), false);
  return sum;
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

// Are any of the bits 1?
bool Bitset::Any() const {
  for (size_t i = 0; i < value_.size(); ++i) {
    if (value_[i]) {
      return true;
    }
  }
  return false;
}

void Bitset::Minorize() {
  assert(value_.size() > 0);
  if (value_[0]) {
    value_.flip();
  }
}

// Copy all of the bits from another bitset into this bitset, starting at
// begin, and optionally flipping the bits as they get copied.
void Bitset::CopyFrom(const Bitset& other, size_t begin, bool flip) {
  assert(begin + other.size() <= size());
  if (flip) {
    for (size_t i = 0; i < other.size(); i++) {
      value_[i + begin] = !other[i];
    }
  } else {
    for (size_t i = 0; i < other.size(); i++) {
      value_[i + begin] = other[i];
    }
  }
}

std::experimental::optional<uint32_t> Bitset::SingletonOption() const {
  bool found_already = false;
  uint32_t found_index;
  for (uint32_t i = 0; i < size(); i++) {
    if (value_[i]) {
      if (found_already) {
        // We previously found an index, so this isn't a singleton.
        return std::experimental::nullopt;
      }  // else
      found_already = true;
      found_index = i;
    }
  }
  if (found_already) {
    return found_index;
  }  // else
  return std::experimental::nullopt;
}

// ** SBN-related functions

// TODO comment
Bitset Bitset::SisterExchange() const {
  assert(size() % 2 == 0);
  Bitset exchanged(size());
  size_t chunk_size = size() / 2;
  for (size_t i = 0; i < chunk_size; i++) {
    exchanged.set(i, value_[i + chunk_size]);
    exchanged.set(i + chunk_size, value_[i]);
  }
  return exchanged;
}

Bitset Bitset::SplitChunk(size_t i) const {
  assert(size() % 2 == 0);
  assert(i < 2);
  size_t chunk_size = size() / 2;
  std::vector<bool> new_value(value_.begin() + int32_t(i * chunk_size),
                              value_.begin() + int32_t((i + 1) * chunk_size));
  return Bitset(new_value);
}

std::string Bitset::ToStringChunked(size_t chunk_count) const {
  assert(size() % chunk_count == 0);
  size_t chunk_size = size() / chunk_count;
  std::string str;
  for (size_t i = 0; i < value_.size(); ++i) {
    str += (value_[i] ? '1' : '0');
    if ((i + 1) % chunk_size == 0 && i + 1 < value_.size()) {
      // The next item will start a new chunk, so add a separator.
      str += '|';
    }
  }
  return str;
}

std::string Bitset::SubsplitToString() const { return ToStringChunked(2); }
std::string Bitset::PCSSToString() const { return ToStringChunked(3); }

size_t Bitset::PCSSChunkSize() const {
  assert(size() % 3 == 0);
  return size() / 3;
}

Bitset Bitset::PCSSChunk(size_t i) const {
  size_t chunk_size = PCSSChunkSize();
  std::vector<bool> new_value(value_.begin() + int32_t(i * chunk_size),
                              value_.begin() + int32_t((i + 1) * chunk_size));
  return Bitset(new_value);
}

Bitset Bitset::PCSSParent() const {
  size_t chunk_size = PCSSChunkSize();
  std::vector<bool> new_value(value_.begin(),
                              value_.begin() + int32_t(2 * chunk_size));
  return Bitset(new_value);
}

Bitset Bitset::PCSSChild() const {
  size_t chunk_size = PCSSChunkSize();
  std::vector<bool> new_value(value_.begin() + int32_t(chunk_size),
                              value_.begin() + int32_t(3 * chunk_size));
  for (size_t i = 0; i < chunk_size; i++) {
    // If A is the child clade, and B is one half of the child split, take the
    // things that are in A but not in B.
    new_value[i] = new_value[i] && !(new_value[i + chunk_size]);
  }
  return Bitset(new_value);
}

bool Bitset::PCSSIsValid() const {
  if (size() % 3 != 0) return false;
  Bitset uncut_parent = PCSSChunk(0);
  Bitset cut_parent = PCSSChunk(1);
  Bitset child = PCSSChunk(2);
  // The parents should be disjoint.
  if ((uncut_parent & cut_parent).Any()) {
    return false;
  }
  // The child should split the cut_parent, so the taxa of child should be a
  // subset of those of cut_parent.
  if ((child & (~cut_parent)).Any()) {
    return false;
  }
  // Something has to be set in each chunk.
  if (!uncut_parent.Any() || !cut_parent.Any() || !child.Any()) {
    return false;
  }
  return true;
}

Bitset Bitset::ChildSubsplit(const Bitset& parent_subsplit,
                             const Bitset& child_half) {
  size_t taxon_count = child_half.size();
  Bitset result(2 * taxon_count);
  assert(result.size() == parent_subsplit.size());
  for (size_t i = 0; i < taxon_count; i++) {
    result.set(i, parent_subsplit[taxon_count + i] ^ child_half[i]);
    result.set(i + taxon_count, child_half[i]);
  }
  return result;
}
