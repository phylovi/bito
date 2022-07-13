// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <cassert>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// TODO: Clean up StrongTypes to match FluentC --
// https://www.fluentcpp.com/2017/05/30/implementing-a-hash-function-for-strong-types/

// // Generic Strong Type (Wrapper for primitive).
// // https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
// template <typename T, typename Parameter>
// class NamedType {
//  public:
//   explicit NamedType(T const &value) : value_(value) {}
//   explicit NamedType(T &&value) : value_(std::move(value)) {}
//   T &get() { return value_; }
//   T const &get() const { return value_; }
//   explicit operator T &() { return value_; }
//   explicit operator T() const { return value_; }

//  private:
//   T value_;
// };

// // Hashable attribute.
// template <typename T>
// struct Hashable {
//   static constexpr bool is_hashable = true;
// };

// // Generic StrongType Hash Function.
// namespace std {
// template <typename T, typename Parameter, typename Converter,
//           template <typename> class... Skills>
// struct hash<NamedTypeImpl<T, Parameter, Converter, Skills...>> {
//   using NamedType = NamedTypeImpl<T, Parameter, Converter, Skills...>;
//   using checkIfHashable = typename std::enable_if<NamedType::is_hashable,
//   void>::type;

//   size_t operator()(NamedTypeImpl<T, Parameter, Converter, Skills...> const &x) const
//   {
//     return std::hash<T>()(x.get());
//   }
// };
// }  // namespace std

// Wrapper for strict typing of primitive types
template <typename T>
struct StrictType {
  explicit StrictType(const T &value) : value_(value) {}

  operator T &() { return value_; }
  operator T() const { return value_; }

  // Outputs Bitset string representation to stream.
  friend std::ostream &operator<<(std::ostream &os, const StrictType<T> &t) {
    os << t.value_;
    return os;
  }

  T value_;
};

struct IdType : public StrictType<size_t> {
  using StrictType<size_t>::StrictType;

  // Outputs string representation to stream.
  friend std::ostream &operator<<(std::ostream &os, const IdType &t) {
    if (t.value_ == NoId) {
      os << "NoId";
    } else {
      os << t.value_;
    }
    return os;
  }

  static constexpr size_t NoId = std::numeric_limits<size_t>::max();
};

// Hash functions for StrictType.
namespace std {
template <>
struct hash<IdType> {
  std::size_t operator()(const IdType &id) const noexcept {
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
  EnumArray(DataType fill_value) { array_.fill(fill_value); }
  EnumArray(std::array<DataType, EnumCount> array) : array_(std::move(array)) {}

  DataType &operator[](const EnumType i) { return array_[static_cast<int>(i)]; }
  const DataType &operator[](const EnumType i) const {
    return array_[static_cast<int>(i)];
  }

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
