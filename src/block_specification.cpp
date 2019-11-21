// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "block_specification.hpp"

void BlockSpecification::Insert(const std::string& key, Coordinates value) {
  SafeInsert(map_, key, value);
}

void BlockSpecification::Insert(const char* key, Coordinates value) {
  Insert(std::string(key), value);
}

BlockSpecification::BlockSpecification(ParamCounts param_counts) {
  size_t next_available_idx = 0;
  for (const auto [block_name, block_size] : param_counts) {
    Insert(block_name, {next_available_idx, block_size});
    next_available_idx += block_size;
  }
  InsertEntireKey({0, next_available_idx});
}

BlockSpecification::Coordinates BlockSpecification::Find(
    const std::string& key) const {
  auto search = map_.find(key);
  if (search == map_.end()) {
    Failwith("Can't find '" + key + "' in block specification!");
  }
  return search->second;
}

void BlockSpecification::InsertEntireKey(Coordinates coordinates) {
  EraseEntireKey();
  Insert(entire_key_, coordinates);
}

void BlockSpecification::Append(const std::string& sub_entire_key,
                                BlockSpecification other) {
  auto next_available_idx = ParameterCount();
  const auto original_next_available_idx = next_available_idx;
  for (const auto [block_name, coordinate] : other.GetMap()) {
    auto [start_idx, block_size] = coordinate;
    if (block_name != entire_key_) {
      auto search = map_.find(block_name);
      if (search != map_.end()) {
        Failwith("Key overlap between BlockSpecifications: " + block_name);
      }  // else
      Insert(block_name, {next_available_idx, block_size});
      next_available_idx += block_size;
    } else {
      Insert(sub_entire_key,
             {start_idx + original_next_available_idx, block_size});
    }
  }
  InsertEntireKey({0, next_available_idx});
}
