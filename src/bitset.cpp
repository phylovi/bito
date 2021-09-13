// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Note that the default constructor for std::vector<bool> is filled with false:
// https://stackoverflow.com/a/22984114/467327

#include "bitset.hpp"

#include <utility>

#include "sugar.hpp"

Bitset::Bitset(std::vector<bool> value) : value_(std::move(value)) {}

Bitset::Bitset(size_t n, bool initial_value) : value_(n, initial_value) {}

Bitset::Bitset(std::string str) : Bitset(str.length()) {
  for (size_t i = 0; i < value_.size(); i++) {
    if (str[i] == '0') {
      value_[i] = false;
    } else if (str[i] == '1') {
      value_[i] = true;
    } else {
      Failwith("String constructor for Bitset must use only 0s or 1s; found '" +
               std::string(1, str[i]) + "'.");
    }
  }
}

bool Bitset::operator[](size_t i) const { return value_[i]; }

size_t Bitset::size() const { return value_.size(); }

void Bitset::set(size_t i, bool value) {
  Assert(i < value_.size(), "i out of range in Bitset::set.");
  value_[i] = value;
}

void Bitset::reset(size_t i) {
  Assert(i < value_.size(), "i out of range in Bitset::reset.");
  value_[i] = false;
}

void Bitset::flip() { value_.flip(); }

bool Bitset::operator==(const Bitset& other) const { return value_ == other.value_; }
bool Bitset::operator!=(const Bitset& other) const { return value_ != other.value_; }
bool Bitset::operator<(const Bitset& other) const { return value_ < other.value_; }
bool Bitset::operator<=(const Bitset& other) const { return value_ <= other.value_; }
bool Bitset::operator>(const Bitset& other) const { return value_ > other.value_; }
bool Bitset::operator>=(const Bitset& other) const { return value_ >= other.value_; }

Bitset Bitset::operator&(const Bitset& other) const {
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator&.");
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] && other.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

Bitset Bitset::operator|(const Bitset& other) const {
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator|.");
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] || other.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

Bitset Bitset::operator^(const Bitset& other) const {
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator^.");
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] != other.value_[i]) {
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
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator&=.");
  for (size_t i = 0; i < value_.size(); i++) {
    value_[i] = value_[i] && other[i];
  }
}

void Bitset::operator|=(const Bitset& other) {
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator|=.");
  for (size_t i = 0; i < value_.size(); i++) {
    value_[i] = value_[i] || other[i];
  }
}

// These methods aren't in the bitset interface.

void Bitset::Zero() { std::fill(value_.begin(), value_.end(), false); }

size_t Bitset::Hash() const { return std::hash<std::vector<bool>>{}(value_); }

std::string Bitset::ToString() const {
  std::string str;
  for (const auto& bit : value_) {
    str += (bit ? '1' : '0');
  }
  return str;
}

bool Bitset::All() const {
  for (const auto& bit : value_) {
    if (!bit) {
      return false;
    }
  }
  return true;
}

bool Bitset::Any() const {
  for (const auto& bit : value_) {
    if (bit) {
      return true;
    }
  }
  return false;
}

bool Bitset::None() const { return !Any(); }

bool Bitset::IsSingleton() const { return SingletonOption().has_value(); }

bool Bitset::IsDisjoint(const Bitset& other) const {
  Assert(size() == other.size(), "Size mismatch in Bitset::IsDisjoint.");
  for (size_t i = 0; i < size(); i++) {
    if (value_[i] && other.value_[i]) {
      return false;
    }
  }
  return true;
}

void Bitset::Minorize() {
  Assert(!(value_.empty()), "Can't Bitset::Minorize an empty bitset.");
  if (value_[0]) {
    value_.flip();
  }
}

void Bitset::CopyFrom(const Bitset& other, size_t begin, bool flip) {
  Assert(begin + other.size() <= size(), "Can't fit copy in Bitset::CopyFrom.");
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

std::optional<uint32_t> Bitset::SingletonOption() const {
  bool found_already = false;
  uint32_t found_index;
  for (uint32_t i = 0; i < size(); i++) {
    if (value_[i]) {
      if (found_already) {
        // We previously found an index, so this isn't a singleton.
        return std::nullopt;
      }  // else
      found_already = true;
      found_index = i;
    }
  }
  if (found_already) {
    return found_index;
  }  // else
  return std::nullopt;
}

size_t Bitset::Count() const { return std::count(value_.begin(), value_.end(), true); }

std::string Bitset::ToIndexSetString() const {
  std::string str;
  for (size_t i = 0; i < size(); i++) {
    if (value_[i]) {
      str += std::to_string(i);
      str += ",";
    }
  }
  if (!str.empty()) {
    str.pop_back();
  }
  return str;
}

// ** SBN-related functions

Bitset Bitset::RotateSubsplit() const {
  Assert(size() % 2 == 0, "Bitset::RotateSubsplit requires an even-size bitset.");
  Bitset exchanged(size());
  size_t chunk_size = size() / 2;
  for (size_t i = 0; i < chunk_size; i++) {
    exchanged.set(i, value_[i + chunk_size]);
    exchanged.set(i + chunk_size, value_[i]);
  }
  return exchanged;
}

std::string Bitset::ToStringChunked(size_t chunk_count) const {
  Assert(size() % chunk_count == 0,
         "Size isn't a multiple of chunk_count in Bitset::ToStringChunked.");
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

std::string Bitset::SubsplitToIndexSetString() const {
  std::string str;
  str += SubsplitChunk(0).ToIndexSetString();
  str += "|";
  str += SubsplitChunk(1).ToIndexSetString();
  return str;
}

size_t Bitset::SubsplitChunkSize() const {
  Assert(size() % 2 == 0, "Size isn't 0 mod 2 in Bitset::SubsplitChunkSize.");
  return size() / 2;
}

Bitset Bitset::SubsplitChunk(size_t i) const {
  Assert(i < 2, "SubsplitChunk index too large.");
  size_t chunk_size = SubsplitChunkSize();
  std::vector<bool> new_value(value_.begin() + int32_t(i * chunk_size),
                              value_.begin() + int32_t((i + 1) * chunk_size));
  return Bitset(std::move(new_value));
}

Bitset Bitset::SubsplitGetChild(size_t i) const {
  Assert(i < 2, "SubsplitGetChild index too large.");
  // The children appear in reverse order; see header file for details.
  Bitset child_1 = SubsplitChunk(0);
  Bitset child_0 = SubsplitChunk(1);
  Assert(child_1 > child_0,
         "Subsplit child 1 should be larger than child 0 in Bitset::SubsplitGetChild.");
  return i == 0 ? child_0 : child_1;
}

bool Bitset::SubsplitIsLeaf() const {
  return SubsplitChunk(0).IsSingleton() && SubsplitChunk(1).None();
}

bool Bitset::SubsplitIsRootsplit() const {
  return SubsplitChunkUnion().All() && SubsplitChunk(0).IsDisjoint(SubsplitChunk(1)) &&
         !SubsplitChunk(0).All();
}

bool Bitset::SubsplitIsRotatedChildOf(const Bitset& other) const {
  return size() == other.size() && SubsplitChunkUnion() == other.SubsplitChunk(0);
}

bool Bitset::SubsplitIsSortedChildOf(const Bitset& other) const {
  return size() == other.size() && SubsplitChunkUnion() == other.SubsplitChunk(1);
}

Bitset Bitset::SubsplitChunkUnion() const {
  Assert(size() % 2 == 0, "Size isn't 0 mod 2 in Bitset::SubsplitChunkUnion.");
  return SubsplitChunk(0) | SubsplitChunk(1);
}

size_t Bitset::PCSPChunkSize() const {
  Assert(size() % 3 == 0, "Size isn't 0 mod 3 in Bitset::PCSPChunkSize.");
  return size() / 3;
}

Bitset Bitset::PCSPChunk(size_t i) const {
  Assert(i < 3, "PCSPChunk index too large.");
  size_t chunk_size = PCSPChunkSize();
  std::vector<bool> new_value(value_.begin() + int32_t(i * chunk_size),
                              value_.begin() + int32_t((i + 1) * chunk_size));
  return Bitset(std::move(new_value));
}

Bitset Bitset::PCSPParentSubsplit() const {
  Bitset sister = PCSPChunk(0);
  Bitset focal = PCSPChunk(1);
  return sister > focal ? sister + focal : focal + sister;
}

Bitset Bitset::PCSPChildSubsplit() const {
  Bitset focal = PCSPChunk(1);
  Bitset child_0 = PCSPChunk(2);
  Bitset child_1(child_0.size());
  for (size_t i = 0; i < child_0.size(); i++) {
    child_1.value_[i] = focal.value_[i] && !child_0.value_[i];
  }
  return child_1 > child_0 ? child_1 + child_0 : child_0 + child_1;
}

std::string Bitset::PCSPToString() const { return ToStringChunked(3); }

bool Bitset::PCSPIsValid() const {
  if (size() % 3 != 0) {
    return false;
  }
  Bitset sister = PCSPChunk(0);
  Bitset focal = PCSPChunk(1);
  Bitset child_0 = PCSPChunk(2);
  // The parents should be disjoint.
  if (!sister.IsDisjoint(focal)) {
    return false;
  }
  // The child should split the focal clade of the parent,
  // so the taxa of child_0 should be a subset of those of focal clade.
  if (!child_0.IsDisjoint(~focal)) {
    return false;
  }
  // Something has to be set in each chunk.
  if (sister.None() || focal.None() || child_0.None()) {
    return false;
  }
  return true;
}

bool Bitset::PCSPIsFake() const {
  Assert(size() % 3 == 0, "Size isn't 0 mod 3 in Bitset::PCSPIsFake.");
  return PCSPChunk(2).None();
}

bool Bitset::PCSPParentIsRootsplit() const {
  Assert(size() % 3 == 0, "Size isn't 0 mod 3 in Bitset::PCSPIsRootsplit.");
  return PCSPParentSubsplit().SubsplitIsRootsplit();
}

SizePair Bitset::PCSPChildSubsplitTaxonCounts() const {
  auto chunk_size = PCSPChunkSize();
  auto total_child_taxon_count =
      std::count(value_.begin() + chunk_size, value_.begin() + 2 * chunk_size, true);
  auto child0_taxon_count =
      std::count(value_.begin() + 2 * chunk_size, value_.end(), true);
  Assert(child0_taxon_count < total_child_taxon_count,
         "PCSPChildSubsplitTaxonCounts: not a proper PCSP bitset.");
  return {static_cast<size_t>(child0_taxon_count),
          static_cast<size_t>(total_child_taxon_count - child0_taxon_count)};
}

Bitset Bitset::Singleton(size_t n, size_t which_on) {
  Assert(which_on < n, "which_on too big in Bitset::Singleton.");
  Bitset singleton(n);
  singleton.set(which_on);
  return singleton;
}

Bitset Bitset::SubsplitOfPair(const Bitset& chunk_0, const Bitset& chunk_1) {
  Assert(chunk_0.IsDisjoint(chunk_1),
         "SubsplitOfPair: given bitsets are not a valid chunk pair.");
  return chunk_1 > chunk_0 ? chunk_1 + chunk_0 : chunk_0 + chunk_1;
}

Bitset Bitset::PCSPOfPair(const Bitset& parent_subsplit, const Bitset& child_subsplit) {
  Assert(
      (child_subsplit.SubsplitIsRotatedChildOf(parent_subsplit) ||
       child_subsplit.SubsplitIsSortedChildOf(parent_subsplit)) &&
          child_subsplit.SubsplitChunk(0).IsDisjoint(child_subsplit.SubsplitChunk(1)),
      "PCSPOfPair: given bitsets are not a valid parent/child pair.");
  Bitset pcsp(3 * parent_subsplit.size() / 2);
  pcsp.CopyFrom((child_subsplit.SubsplitIsRotatedChildOf(parent_subsplit)
                     ? parent_subsplit.RotateSubsplit()
                     : parent_subsplit),
                0, false);
  pcsp.CopyFrom(child_subsplit.SubsplitGetChild(0), parent_subsplit.size(), false);
  return pcsp;
}

Bitset Bitset::FakeSubsplit(const Bitset& nonzero_contents) {
  Bitset fake(2 * nonzero_contents.size());
  // Put the nonzero contents on the left of the fake subsplit.
  fake.CopyFrom(nonzero_contents, 0, false);
  return fake;
}

void AssertSisterAndLeafSubsplit(const Bitset& subsplit) {
  Assert(subsplit.SubsplitChunk(0).Any() && subsplit.SubsplitChunk(1).IsSingleton(),
         "Assertion failed: we want the left-hand chunk of the subsplit be "
         "non-empty and the right-hand chunk be a singleton.");
}

Bitset Bitset::FakeChildSubsplit(const Bitset& parent_subsplit) {
  AssertSisterAndLeafSubsplit(parent_subsplit);
  // Put the right-hand chunk of the subsplit as the nonzero contents of the fake
  // subsplit.
  return FakeSubsplit(parent_subsplit.SubsplitChunk(1));
}

Bitset Bitset::FakePCSP(const Bitset& parent_subsplit) {
  AssertSisterAndLeafSubsplit(parent_subsplit);
  const auto taxon_count = parent_subsplit.SubsplitChunkSize();
  Bitset fake(3 * taxon_count);
  // Put the nonzero contents on the left of the fake subsplit.
  fake.CopyFrom(parent_subsplit, 0, false);
  return fake;
}

Bitset Bitset::DAGRootSubsplitOfTaxonCount(const size_t taxon_count) {
  Bitset zeros(taxon_count);
  return ~zeros + zeros;
}

Bitset Bitset::RootsplitOfHalf(const Bitset& subsplit_half) {
  Bitset half = subsplit_half;
  half.Minorize();
  return ~half + half;
}

Bitset Bitset::PCSPOfRootsplit(const Bitset& rootsplit) {
  Assert(rootsplit.SubsplitIsRootsplit(),
         "Given subsplit is not rootsplit in Bitset::PCSPOfRootsplit.");
  return PCSPOfPair(DAGRootSubsplitOfTaxonCount(rootsplit.size() / 2), rootsplit);
}

Bitset Remap(const Bitset& bitset, const SizeOptionVector& idx_table) {
  Bitset result(idx_table.size(), false);
  for (size_t i = 0; i < idx_table.size(); ++i) {
    if (idx_table[i].has_value() && bitset[idx_table[i].value()]) {
      result.set(i, true);
    }
  }
  return result;
}
