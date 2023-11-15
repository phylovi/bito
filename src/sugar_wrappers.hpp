// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Templates for a strong typing. Creates a lightweight wrapper for primitive types.
// Enforces explicit casting between differenct strong types and the underlying
// primitive. Also includes template for inheriting primitive's hashing function for use
// in stl types.
// https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
//
// Templates for enumerated types. Creates a wrapper class for a base enum class.
// Includes using enums to access stl storage classes.

#pragma once

#include <cassert>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <type_traits>
#include <iterator>

// ** Strongly Typed Wrappers

// Wrapper for ID types.
template <typename TypeNameTag>
struct GenericId {
  using SelfType = GenericId<TypeNameTag>;
  using UnderlyingType = size_t;
  static constexpr size_t NoId = std::numeric_limits<size_t>::max();

  GenericId() = default;
  explicit GenericId(UnderlyingType const &value) : value_(value) {}
  template <typename T_ = UnderlyingType>
  GenericId(
      UnderlyingType &&value,
      typename std::enable_if<!std::is_reference<T_>{}, std::nullptr_t>::type = nullptr)
      : value_(std::move(value)) {}

  // Explicit cast to primitive.
  explicit operator UnderlyingType &() { return value_; }
  explicit operator UnderlyingType() const { return value_; }

  // Compare to its own type.
  int Compare(const SelfType &other) const { return Compare(other.value_); }
  bool operator==(const SelfType &other) const { return value_ == other.value_; }
  bool operator!=(const SelfType &other) const { return value_ != other.value_; }
  bool operator>(const SelfType &other) const { return value_ > other.value_; }
  bool operator<(const SelfType &other) const { return value_ < other.value_; }
  bool operator>=(const SelfType &other) const { return value_ >= other.value_; }
  bool operator<=(const SelfType &other) const { return value_ <= other.value_; }
  // Compare to its primitive.
  int Compare(const UnderlyingType &other) const {
    return (value_ > other) ? value_ - other : other - value_;
  }
  bool operator==(const UnderlyingType &other) const { return value_ == other; }
  bool operator!=(const UnderlyingType &other) const { return value_ != other; }
  bool operator>(const UnderlyingType &other) const { return value_ > other; }
  bool operator<(const UnderlyingType &other) const { return value_ < other; }
  bool operator>=(const UnderlyingType &other) const { return value_ >= other; }
  bool operator<=(const UnderlyingType &other) const { return value_ <= other; }

  // Increment
  SelfType &operator++() {
    value_++;
    return *this;
  }
  SelfType operator++(int) {
    SelfType temp = *this;
    ++*this;
    return temp;
  }

  // I/O
  static std::string PrefixToString() { return "Id"; }
  std::string ToString(const bool include_prefix = true) const {
    std::stringstream os;
    os << (include_prefix ? PrefixToString() : "")
       << "::" << ((value_ == NoId) ? "NoId" : std::to_string(value_));
    return os.str();
  }
  friend std::ostream &operator<<(std::ostream &os, const SelfType &obj) {
    os << obj.ToString(true);
    return os;
  }

  UnderlyingType value_ = NoId;
};

constexpr size_t NoId = std::numeric_limits<size_t>::max();

// Generic hash function for GenericIds.
namespace std {
template <typename TypeNameTag>
struct hash<GenericId<TypeNameTag>> {
  std::size_t operator()(const GenericId<TypeNameTag> &id) const noexcept {
    std::size_t value_hash = std::hash<size_t>()(id.value_);
    return value_hash;
  }
};
}  // namespace std

// ** Enumerated Type Wrappers

// Generic Iterator for enum class types.
// Requires that enum's underlying types have no gaps.
template <typename EnumType, EnumType FirstEnum, EnumType LastEnum>
class EnumIterator {
  typedef typename std::underlying_type<EnumType>::type val_t;
  int val;

 public:
  EnumIterator(const EnumType &f) : val(static_cast<val_t>(f)) {}
  EnumIterator() : val(static_cast<val_t>(FirstEnum)) {}
  EnumIterator operator++() {
    ++val;
    return *this;
  }
  EnumType operator*() { return static_cast<EnumType>(val); }
  EnumIterator begin() { return *this; }  // default ctor is good
  EnumIterator end() {
    static const EnumIterator endIter = ++EnumIterator(LastEnum);  // cache it
    return endIter;
  }
  bool operator!=(const EnumIterator &i) { return val != i.val; }
};

// Generic Array for using class enum for index access.
template <class EnumType, size_t EnumCount, class DataType>
class EnumArray {
 public:
  EnumArray() = default;
  EnumArray(std::array<DataType, EnumCount> array) : array_(std::move(array)) {}

  DataType &operator[](const EnumType i) { return array_[static_cast<int>(i)]; }
  const DataType &operator[](const EnumType i) const {
    return array_[static_cast<int>(i)];
  }

  void fill(DataType fill_value) { array_.fill(fill_value); }
  int size() const { return array_.size(); }

  auto begin() const { return std::begin(array_); }
  auto end() const { return std::end(array_); }

  friend std::ostream &operator<<(std::ostream &os, const EnumArray &array) {
    os << array.array_;
    return os;
  }

 private:
  std::array<DataType, EnumCount> array_;
};

// Generic Wrapper with collection of static functions for using class enum with
// common data structures.
template <class EnumType, class UnderlyingType, size_t EnumCount, EnumType FirstEnum,
          EnumType LastEnum>
class EnumWrapper {
 public:
  using Type = EnumType;
  using Iterator = EnumIterator<EnumType, FirstEnum, LastEnum>;
  template <class DataType>
  using Array = EnumArray<EnumType, EnumCount, DataType>;

  static inline const EnumType First = FirstEnum;
  static inline const EnumType Last = LastEnum;
  static inline const size_t Count = EnumCount;

  static UnderlyingType GetIndex(const EnumType i) {
    return static_cast<UnderlyingType>(i);
  }

  static std::array<Type, Count> TypeArray() {
    std::array<Type, Count> arr;
    size_t i = 0;
    for (const auto e : Iterator()) {
      arr[i++] = e;
    }
    return arr;
  }

  static std::string ToString(const EnumType i) {
    std::stringstream os;
    os << "Enum::" << std::to_string(GetIndex(i));
    return os.str();
  }
};
