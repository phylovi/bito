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

#include "intpack.hpp"
#include "prettyprint.hpp"

// Put typedefs that are built of STL types here.
using Tag = uint64_t;
using SymbolVector = std::vector<int>;
using SizeVector = std::vector<size_t>;
using SizeVectorVector = std::vector<SizeVector>;
using DoubleVectorVector = std::vector<std::vector<double>>;
using TagDoubleMap = std::unordered_map<Tag, double>;
using TagSizeMap = std::unordered_map<Tag, size_t>;
using TagStringMap = std::unordered_map<Tag, std::string>;
using StringStringMap = std::unordered_map<std::string, std::string>;
using CharIntMap = std::unordered_map<char, int>;
using StringSizeMap = std::unordered_map<std::string, size_t>;
using StringDoubleMap = std::unordered_map<std::string, double>;
using DoubleVectorOption = std::optional<std::vector<double>>;
using TagStringMapOption = std::optional<TagStringMap>;
using StringVector = std::vector<std::string>;
using StringVectorVector = std::vector<StringVector>;
using StringSet = std::unordered_set<std::string>;
using StringSetVector = std::vector<StringSet>;
using StringDoubleVector = std::vector<std::pair<std::string, double>>;
using SizeDoubleVectorMap = std::unordered_map<size_t, std::vector<double>>;
using DoublePair = std::pair<double, double>;
using SizePair = std::pair<size_t, size_t>;
using SizeOptionVector = std::vector<std::optional<size_t>>;

inline uint32_t MaxLeafIDOfTag(Tag tag) { return UnpackFirstInt(tag); }
inline uint32_t LeafCountOfTag(Tag tag) { return UnpackSecondInt(tag); }

// We implement problems in terms of exceptions. That means that they work great
// in Jupyter notebooks.
//
// This macro always evaluates the argument. We use a macro for the stupid
// reason that then the assert can go away upon using NDEBUG.
#ifdef NDEBUG
#define Assert(to_evaluate, message) ((void)(to_evaluate));
#else
#define Assert(to_evaluate, message) \
  ({                                 \
    if (!(to_evaluate)) {            \
      Failwith(message);             \
    }                                \
  })
#endif
// Use Failwith when it's a problem with input data versus a problem
// with program logic. That way we can turn off Assert if we want to.
// As you can see Assert is implemented in terms of Failwith.
//
// Here we use a macro to avoid "control may reach end of non-void function"
// errors. We shouldn't have to return when we throw an exception.
#define Failwith(message)                         \
  ({                                              \
    std::string str_message(message);             \
    str_message.append(" (");                     \
    str_message.append(__FILE__);                 \
    str_message.append(":");                      \
    str_message.append(std::to_string(__LINE__)); \
    str_message.append(" in ");                   \
    str_message.append(__func__);                 \
    str_message.append(")");                      \
    throw std::runtime_error(str_message);        \
  })

template <class Key, class T, class Hash>
constexpr void SafeInsert(std::unordered_map<Key, T, Hash> &map, const Key &k,
                          const T &v) {
  Assert(map.insert({k, v}).second, "Failed map insertion!");
}

template <class Key, class T, class Hash>
constexpr void SafeInsert(std::unordered_map<Key, T, Hash> &map, Key &&k, T &&v) {
  Assert(map.insert({k, v}).second, "Failed map insertion!");
}

template <class Key, class T, class Hash>
constexpr void SafeInsert(std::map<Key, T, Hash> &map, const Key &k, const T &v) {
  Assert(map.insert({k, v}).second, "Failed map insertion!");
}

template <class Key, class T, class Hash>
constexpr void SafeInsert(std::map<Key, T, Hash> &map, Key &&k, T &&v) {
  Assert(map.insert({k, v}).second, "Failed map insertion!");
}

template <class Key, class Hash>
constexpr void SafeInsert(std::unordered_set<Key, Hash> &set, const Key &k) {
  Assert(set.insert(k).second, "Failed set insertion!");
}

template <class Key, class Hash>
constexpr void SafeInsert(std::unordered_set<Key, Hash> &set, Key &&k) {
  Assert(set.insert(k).second, "Failed set insertion!");
}

template <class Key, class T, class Hash>
T AtWithDefault(const std::unordered_map<Key, T, Hash> &map, const Key &key,
                T default_value) {
  auto search = map.find(key);
  if (search == map.end()) {
    return default_value;
  }
  return search->second;
}

template <typename T>
std::string GenericToString(T x) {
  std::stringstream ss;
  ss << x << std::endl;
  return ss.str();
}

template <class Key, class T>
std::unordered_map<Key, T> UnorderedMapOf(const std::vector<std::pair<Key, T>> &v) {
  std::unordered_map<Key, T> m;
  for (const auto &[key, value] : v) {
    SafeInsert(m, key, value);
  }
  return m;
}
