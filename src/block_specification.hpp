// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// This class describes the structure of parameter collections that sit in
// contiguous blocks. These structures are maps from string keys (which are the
// block names) to Coordinates. See the unit tests at the bottom to see how they
// work.
//
// There is a special key called entire_key_ which gives the entire span of the
// block.

#ifndef SRC_BLOCK_SPECIFICATION_HPP_
#define SRC_BLOCK_SPECIFICATION_HPP_

#include "sugar.hpp"

class BlockSpecification {
 public:
  // Block Coordinates are (starting index, block size) in a tuple.
  using Coordinates = std::tuple<size_t, size_t>;
  using ParamCounts = std::vector<std::tuple<std::string, size_t>>;
  using UnderlyingMapType = std::unordered_map<std::string, Coordinates>;

  // Given a map of block names to the number of parameters they have, here we
  // build out a block specification.
  BlockSpecification(ParamCounts param_counts);
  const UnderlyingMapType& GetMap() const { return map_; }

  Coordinates Find(const std::string& key) const;

  void Insert(const std::string& key, Coordinates value) {
    SafeInsert(map_, key, value);
  }

  // Insert for string literals
  void Insert(const char* key, Coordinates value) {
    Insert(std::string(key), value);
  }

  // See top for description of what Entire means in this case.
  void InsertEntireKey(Coordinates coordinates);
  void EraseEntireKey() { map_.erase(entire_key_); }

  // Incorporate one BlockSpecification into `this` one, starting at
  // next_available_index which we mutate along the way. The "entire" block
  // coordinates from other get incorporated into this with key sub_entire_key,
  // with the coordinates incremented by the next_available_idx as passed into
  // this function.
  void Append(const std::string& sub_entire_key, BlockSpecification other);

  // The complete range of parameter counts.
  size_t ParameterCount() const { return std::get<1>(Find(entire_key_)); };

  inline const static std::string entire_key_ = "entire";

 private:
  UnderlyingMapType map_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("BlockSpecification") {
  // As an example, kazoo has 4 parameters, and jordan has 23.
  BlockSpecification spec({{"kazoo", 4}, {"jordan", 23}});
  // Here we can see that the specification stores the starting index and then
  // the number of parameters.
  CHECK_EQ(spec.Find("kazoo"), BlockSpecification::Coordinates({0, 4}));
  CHECK_EQ(spec.Find("jordan"), BlockSpecification::Coordinates({4, 23}));
  spec.Append("entire turbo boost",
              BlockSpecification({{"turbo", 666}, {"boost", 42}}));
  // After appending, we find boost at 23+4=27.
  CHECK_EQ(spec.Find("turbo"), BlockSpecification::Coordinates({27, 666}));
  // Turbo is at 666+27 = 693.
  CHECK_EQ(spec.Find("boost"), BlockSpecification::Coordinates({693, 42}));
  // The entire turbo boost starts at 27 and is 42+666 = 708 long.
  CHECK_EQ(spec.Find("entire turbo boost"),
           BlockSpecification::Coordinates({27, 708}));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BLOCK_SPECIFICATION_HPP_
