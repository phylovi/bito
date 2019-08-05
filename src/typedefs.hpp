// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TYPEDEFS_HPP_
#define SRC_TYPEDEFS_HPP_

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "optional.hpp"

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

#endif  // SRC_TYPEDEFS_HPP_
