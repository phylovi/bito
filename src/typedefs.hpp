// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TYPEDEFS_HPP_
#define SRC_TYPEDEFS_HPP_

#include <optional>
#include <unordered_map>

// Put typedefs that are built of STL types here.
typedef uint64_t Tag;
typedef std::vector<int> SymbolVector;
typedef std::unordered_map<Tag, double> TagDoubleMap;
typedef std::optional<TagDoubleMap> TagDoubleMapOption;
typedef std::unordered_map<Tag, std::string> TagStringMap;
typedef std::optional<TagStringMap> TagStringMapOption;
typedef std::unordered_map<std::string, std::string> StringStringMap;
typedef std::unordered_map<char, int> CharIntMap;

#endif  // SRC_TYPEDEFS_HPP_
