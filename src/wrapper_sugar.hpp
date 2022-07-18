// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Templates for a strong typing. Creates a lightweight wrapper for primitive types.
// Enforces explicit casting between differenct strong types and the underlying
// primitive. Also includes template for inheriting primitive's hashing function for use
// in stl types. For more information, read:
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

// Is Comparable

// Wrapper for strong typing of primitive types
template <typename Type, typename TypeNameTag>
struct StrongType {
  using UnderlyingType = Type;

  explicit StrongType(UnderlyingType const &value) : value_(value) {}
  // Template constructor for reference type wrappers.
  template <typename T_ = UnderlyingType>
  explicit StrongType(
      UnderlyingType &&value,
      typename std::enable_if<!std::is_reference<T_>{}, std::nullptr_t>::type = nullptr)
      : value_(std::move(value)) {}

  explicit StrongType(UnderlyingType &&value) : value_(std::move(value)) {}

  UnderlyingType &get() { return value_; }
  UnderlyingType const &get() const { return value_; }
  operator UnderlyingType &() { return value_; }
  operator UnderlyingType() const { return value_; }

  // Can implicitly assign GenericId to size_t.
  StrongType<Type, TypeNameTag> &operator=(const UnderlyingType &new_value) {
    value_ = new_value;
    return *this;
  }

  // Compare to its own type.
  int Compare(const StrongType<Type, TypeNameTag> &other) const {
    return Compare(other.value_);
  }
  bool operator==(const StrongType<Type, TypeNameTag> &other) const {
    return Compare(other) == 0;
  }
  bool operator>(const StrongType<Type, TypeNameTag> &other) const {
    return Compare(other) > 0;
  }
  bool operator<(const StrongType<Type, TypeNameTag> &other) const {
    return Compare(other) < 0;
  }
  // Compare to its underlying type.
  int Compare(const UnderlyingType &other) const { return value_ - other; }
  bool operator==(const UnderlyingType &other) const { return Compare(other) == 0; }
  bool operator>(const UnderlyingType &other) const { return Compare(other) > 0; }
  bool operator<(const UnderlyingType &other) const { return Compare(other) < 0; }

  // Outputs string representation to stream.
  friend std::ostream &operator<<(std::ostream &os,
                                  const StrongType<Type, TypeNameTag> &obj) {
    os << obj.value_;
    return os;
  }

  Type value_;
};

// Generic hash function for StrongType.
namespace std {
template <typename Type, typename TypeNameTag>
struct hash<StrongType<Type, TypeNameTag>> {
  std::size_t operator()(const StrongType<Type, TypeNameTag> &id) const noexcept {
    std::size_t value_hash = std::hash<Type>()(id.value_);
    return value_hash;
  }
};
}  // namespace std

//
template <typename TypeNameTag>
struct GenericId : public StrongType<size_t, TypeNameTag> {
  using NameTag = TypeNameTag;
  using UnderlyingType = size_t;

  explicit GenericId() : StrongType(NoId) {}
  explicit GenericId(UnderlyingType const &value) : StrongType(value) {}
  // Template constructor for reference type wrappers.
  template <typename T_ = UnderlyingType>
  explicit GenericId(
      UnderlyingType &&value,
      typename std::enable_if<!std::is_reference<T_>{}, std::nullptr_t>::type = nullptr)
      : StrongType(std::move(value)) {}

  static std::string PrefixToString() { return "Id"; }

  std::string ToString(const bool include_prefix = true) const {
    std::stringstream os;
    os << (include_prefix ? PrefixToString() : "")
       << "::" << ((value_ == NoId) ? "NoId" : std::to_string(value_));
    return os.str();
  }

  // Can implicitly assign GenericId to size_t.
  GenericId<TypeNameTag> &operator=(const UnderlyingType &new_value) {
    value_ = new_value;
    return *this;
  }

  // Compare to its own type.
  int Compare(const GenericId<TypeNameTag> &other) const {
    return Compare(other.value_);
  }
  bool operator==(const GenericId<TypeNameTag> &other) const {
    return Compare(other) == 0;
  }
  bool operator>(const GenericId<TypeNameTag> &other) const {
    return Compare(other) > 0;
  }
  bool operator<(const GenericId<TypeNameTag> &other) const {
    return Compare(other) < 0;
  }

  friend std::ostream &operator<<(std::ostream &os, const GenericId<TypeNameTag> &obj) {
    os << obj.ToString(true);
    return os;
  }

  static constexpr size_t NoId = std::numeric_limits<size_t>::max();
};

constexpr size_t NoId = GenericId<struct NoNameTag>::NoId;

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
