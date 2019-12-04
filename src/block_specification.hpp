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

#include "eigen_sugar.hpp"
#include "sugar.hpp"

class BlockSpecification {
 public:
  // Block Coordinates are (starting index, block size) in a tuple.
  using Coordinates = std::tuple<size_t, size_t>;
  // ParamCounts are vectors of tuples of (block name, number of parameters for
  // that block).
  using ParamCounts = std::vector<std::tuple<std::string, size_t>>;
  using UnderlyingMapType = std::unordered_map<std::string, Coordinates>;

  // These are handy structures: maps from the block specification keys to
  // segments (i.e. sub-vectors) and blocks (i.e. sub-matrices). We can then
  // read and write to these values, which will be reflected in the original
  // parameter vector/matrix.
  using ParameterSegmentMap = std::unordered_map<std::string, EigenVectorXdRef>;
  using ParameterBlockMap = std::unordered_map<std::string, EigenMatrixXdRef>;

  BlockSpecification(ParamCounts param_counts);

  Coordinates Find(const std::string& key) const;
  void Insert(const std::string& key, Coordinates value);
  // Insert for string literals.
  void Insert(const char* key, Coordinates value);

  // Incorporate one BlockSpecification into `this` one. The "entire" block
  // coordinates from other get incorporated into this with key sub_entire_key.
  // The "entire" block coordinates are then updated.
  void Append(const std::string& sub_entire_key, BlockSpecification other);

  void CheckParameterVectorSize(const EigenVectorXdRef param_vector) const;
  void CheckParameterMatrixSize(const EigenMatrixXdRef param_matrix) const;

  // These methods allow us to pull out segments (i.e. sub-vectors) from vectors
  // and blocks from matrices depending on the coordinates of the block
  // specification. They are very useful for writing the SetParameters method of
  // BlockModels.
  EigenVectorXdRef ExtractSegment(EigenVectorXdRef param_vector,
                                  std::string key) const;
  EigenMatrixXdRef ExtractBlock(EigenMatrixXdRef param_matrix,
                                std::string key) const;

  // These are explained in the definition of ParameterSegmentMap and
  // ParameterBlockMap.
  ParameterSegmentMap ParameterSegmentMapOf(
      EigenVectorXdRef param_vector) const;
  ParameterBlockMap ParameterBlockMapOf(EigenMatrixXdRef param_matrix) const;

  const UnderlyingMapType& GetMap() const { return map_; }
  // The complete range of parameter counts.
  size_t ParameterCount() const { return std::get<1>(Find(entire_key_)); };

  // Here's how we set the "entire" key. In C++ we can just use these key
  // variables, but in Python we use the strings as dictionary keys.
  inline const static std::string entire_key_ = "entire";

 private:
  UnderlyingMapType map_;

  // See top for description of what Entire means in this case.
  void InsertEntireKey(Coordinates coordinates);
  void EraseEntireKey() { map_.erase(entire_key_); }
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
  auto [turbo_index, turbo_size] = spec.Find("turbo");
  // After appending, we find turbo at 23+4=27 or 23+4+666 = 693 depending on
  // the order of the unordered_map.
  CHECK((turbo_index == 27 || turbo_index == 693));
  CHECK_EQ(turbo_size, 666);
  auto [boost_index, boost_size] = spec.Find("boost");
  CHECK((boost_index == 27 || boost_index == 693));
  CHECK_EQ(boost_size, 42);
  // The entire turbo boost starts at 27 and is 42+666 = 708 long.
  CHECK_EQ(spec.Find("entire turbo boost"),
           BlockSpecification::Coordinates({27, 708}));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BLOCK_SPECIFICATION_HPP_
