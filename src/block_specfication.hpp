// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BLOCK_SPECIFICATION_HPP_
#define SRC_BLOCK_SPECIFICATION_HPP_

#include "sugar.hpp"

using BlockCoordinates = std::tuple<size_t, size_t>;

class BlockSpecification {
 public:
  using UnderlyingMapType = std::unordered_map<std::string, BlockCoordinates>;

  // This allows us to use all of the initializer lists for unordered_map to
  // make BlockSpecifications.
  BlockSpecification(UnderlyingMapType init) : map_(std::move(init)) {}

  const UnderlyingMapType& GetMap() const { return map_; }

  BlockCoordinates Find(const std::string& key) {
    auto search = map_.find(key);
    if (search == map_.end()) {
      Failwith("Can't find '" + key + "' in block specification!");
    }
    return search->second;
  }

  void Insert(const std::string& key, BlockCoordinates value) {
    SafeInsert(map_, key, value);
  }

  // Insert for string literals
  void Insert(const char* key, BlockCoordinates value) {
    Insert(std::string(key), value);
  }

  // The complete range of parameter counts.
  size_t ParameterCount() { return std::get<1>(Find(entire_key_)); };

  inline const static std::string entire_key_ = "entire";

 private:
  UnderlyingMapType map_;
};

#endif  // SRC_BLOCK_SPECIFICATION_HPP_
