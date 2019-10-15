// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SUGAR_HPP_
#define SRC_SUGAR_HPP_

#include <cassert>
#include <optional>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "prettyprint.hpp"

// Put typedefs that are built of STL types here.
typedef uint64_t Tag;
typedef std::vector<int> SymbolVector;
typedef std::vector<size_t> SizeVector;
typedef std::vector<SizeVector> SizeVectorVector;
using DoubleVectorVector = std::vector<std::vector<double>>;
typedef std::unordered_map<Tag, double> TagDoubleMap;
typedef std::unordered_map<Tag, size_t> TagSizeMap;
typedef std::unordered_map<Tag, std::string> TagStringMap;
typedef std::unordered_map<std::string, std::string> StringStringMap;
typedef std::unordered_map<char, int> CharIntMap;
// This will be STL in C++17 but we include the above header to fake it.
typedef std::optional<std::vector<double>> DoubleVectorOption;
typedef std::optional<TagStringMap> TagStringMapOption;
typedef std::vector<std::string> StringVector;
typedef std::vector<StringVector> StringVectorVector;
typedef std::unordered_set<std::string> StringSet;
typedef std::vector<StringSet> StringSetVector;
typedef std::unordered_map<std::string, size_t> StringSizeMap;

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

template <class Key, class Hash>
constexpr void SafeInsert(std::unordered_set<Key, Hash> &set, const Key &k) {
  Assert(set.insert(k).second, "Failed set insertion!");
}

#endif  // SRC_SUGAR_HPP_
