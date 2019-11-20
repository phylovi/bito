// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BLOCK_SPECIFICATION_HPP_
#define SRC_BLOCK_SPECIFICATION_HPP_

#include "sugar.hpp"

class BlockSpecification {
 public:
  // The coordinates of a block consist of the starting index and the block
  // size.
  using Coordinates = std::tuple<size_t, size_t>;
  using ParamCounts = std::unordered_map<std::string, size_t>;
  using UnderlyingMapType = std::unordered_map<std::string, Coordinates>;

  // Given a map of block names to the number of parameters they have, here we
  // build out a block specification.
  BlockSpecification(ParamCounts param_counts) {
    size_t next_available_idx = 0;
    for (const auto [block_name, block_size] : param_counts) {
      Insert(block_name, {next_available_idx, block_size});
      next_available_idx += block_size;
    }
    InsertEntireKey({0, next_available_idx});
  }

  const UnderlyingMapType& GetMap() const { return map_; }

  Coordinates Find(const std::string& key) const {
    auto search = map_.find(key);
    if (search == map_.end()) {
      Failwith("Can't find '" + key + "' in block specification!");
    }
    return search->second;
  }

  void Insert(const std::string& key, Coordinates value) {
    SafeInsert(map_, key, value);
  }

  // Insert for string literals
  void Insert(const char* key, Coordinates value) {
    Insert(std::string(key), value);
  }

  void InsertEntireKey(Coordinates coordinates) {
    EraseEntireKey();
    Insert(entire_key_, coordinates);
  }
  void EraseEntireKey() { map_.erase(entire_key_); }

  // Incorporate one BlockSpecification into this, starting at
  // next_available_index which we mutate along the way. The "entire" block
  // coordinates from other get incorporated into this with key sub_entire_key,
  // with the coordinates incremented by the next_available_idx as passed into
  // this function.
  void Append(const std::string& sub_entire_key, BlockSpecification other) {
    auto next_available_idx = ParameterCount();
    const auto original_next_available_idx = next_available_idx;
    for (const auto [block_name, coordinate] : other.GetMap()) {
      auto [start_idx, block_size] = coordinate;
      if (block_name != entire_key_) {
        Insert(block_name, {next_available_idx, block_size});
        next_available_idx += block_size;
      } else {
        Insert(sub_entire_key,
               {start_idx + original_next_available_idx, block_size});
      }
    }
    InsertEntireKey({0, next_available_idx});
  }

  // The complete range of parameter counts.
  size_t ParameterCount() const { return std::get<1>(Find(entire_key_)); };

  inline const static std::string entire_key_ = "entire";

 private:
  UnderlyingMapType map_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("BlockSpecification") {
  // As an example, kazoo has 4 parameters, and jordan has 23.
  BlockSpecification spec({{"jordan", 23}, {"kazoo", 4}});
  // Here we can see that the specification stores the starting index and then
  // the number of parameters.
  // Note that the order of the blocks is laid out by the (arbitrary) order of
  // keys in unordered_map, not in terms of their order in the initializer list.
  CHECK_EQ(spec.Find("kazoo"), BlockSpecification::Coordinates({0, 4}));
  CHECK_EQ(spec.Find("jordan"), BlockSpecification::Coordinates({4, 23}));
  spec.Append("entire turbo boost",
              BlockSpecification({{"turbo", 42}, {"boost", 666}}));
  // After appending, we find boost at 23+4=27.
  CHECK_EQ(spec.Find("boost"), BlockSpecification::Coordinates({27, 666}));
  // Turbo is at 666+27 = 693.
  CHECK_EQ(spec.Find("turbo"), BlockSpecification::Coordinates({693, 42}));
  // The entire turbo boost starts at 27 and is 42+666 = 708 long.
  CHECK_EQ(spec.Find("entire turbo boost"),
           BlockSpecification::Coordinates({27, 708}));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BLOCK_SPECIFICATION_HPP_
