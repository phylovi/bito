// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SUGAR_HPP_
#define SRC_SUGAR_HPP_

#include <cassert>
#include <experimental/optional>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Put typedefs that are built of STL types here.
typedef uint64_t Tag;
typedef std::vector<int> SymbolVector;
typedef std::vector<size_t> SizeVector;
typedef std::vector<SizeVector> SizeVectorVector;
typedef std::unordered_map<Tag, double> TagDoubleMap;
typedef std::unordered_map<Tag, size_t> TagSizeMap;
typedef std::unordered_map<Tag, std::string> TagStringMap;
typedef std::unordered_map<std::string, std::string> StringStringMap;
typedef std::unordered_map<char, int> CharIntMap;
// This will be STL in C++17 but we include the above header to fake it.
typedef std::experimental::optional<std::vector<double>> DoubleVectorOption;
typedef std::experimental::optional<TagStringMap> TagStringMapOption;
typedef std::vector<std::string> StringVector;
typedef std::unordered_set<std::string> StringSet;
typedef std::vector<StringSet> StringSetVector;

// Here we use bad style by defining a macro. This is in preparation for when we
// can use something nicer in a future C++.
#define Assert(status, message) assert(status &&message)
inline void Failwith(const std::string &message) {
  std::cerr << message << std::endl;
  abort();
}

template <class Key, class T, class Hash>
constexpr void SafeInsert(std::unordered_map<Key, T, Hash> &map, const Key &k,
                          const T &v) {
  auto status = map.insert({k, v}).second;
  Assert(status, "Failed map insertion!");
}

template <class Key, class Hash>
constexpr void SafeInsert(std::unordered_set<Key, Hash> &set, const Key &k) {
  auto status = set.insert(k).second;
  Assert(status, "Failed set insertion!");
}

#endif  // SRC_SUGAR_HPP_
