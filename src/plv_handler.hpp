// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// PLVHandler is used by GPOperations to get index location of PLVs stored on the
// GPEngine.

#pragma once

#include "sugar.hpp"

class PLVHandler {
 public:
  // We store 6 PLVs per subsplit, and index them according to this enum. (See
  // GPEngine::plvs_).
  static inline const size_t plv_count_ = 6;
  enum class PLVType : size_t {
    P,          // p(s)
    PHatRight,  // phat(s_right)
    PHatLeft,   // phat(s_left)
    RHat,       // rhat(s_right) = rhat(s_left)
    RRight,     // r(s_right)
    RLeft,      // r(s_left)
  };
  typedef EnumIterator<PLVType, PLVType::P, PLVType::RLeft> PLVTypeIterator;

  // Get the `GPEngine::plvs_` index of given node's given PLV type from DAG with
  // node_count.
  static size_t GetPLVIndex(const PLVType plv_type, const size_t node_count,
                            const size_t node_idx) {
    return GetPLVIndex(GetPLVTypeIndex(plv_type), node_count, node_idx);
  };

  static PLVType RPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::RLeft : PLVType::RRight;
  };

  static PLVType PPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::PHatLeft : PLVType::PHatRight;
  };

  // ** I/O

  static std::string PLVTypeToString(const PLVType plv_type) {
    return "PLVType::" + std::to_string(GetPLVTypeIndex(plv_type));
  };

  friend std::ostream &operator<<(std::ostream &os, const PLVType plv_type) {
    os << PLVTypeToString(plv_type);
    return os;
  };

 protected:
  static size_t GetPLVTypeIndex(const PLVType plv_type) {
    return static_cast<std::underlying_type<PLVType>::type>(plv_type);
  };

  static size_t GetPLVIndex(const size_t plv_type_idx, const size_t node_count,
                            const size_t node_idx) {
    return (plv_type_idx * node_count) + node_idx;
  };
};

#ifdef DOCTEST_LIBRARY_INCLUDED

using PLVType = PLVHandler::PLVType;

// Check that PLV iterator iterates over all PLVs exactly once.
TEST_CASE("PLVHandler: EnumIterator") {
  const std::array<PLVType, PLVHandler::plv_count_> plv_types{
      PLVType::P,    PLVType::PHatRight, PLVType::PHatLeft,
      PLVType::RHat, PLVType::RRight,    PLVType::RLeft};
  std::map<PLVType, size_t> plv_visited_map;
  // Iterate using vector.
  for (const PLVType plv_type : plv_types) {
    plv_visited_map.insert({plv_type, 0});
  }
  // Iterate using EnumIterator.
  for (const PLVType plv_type : PLVHandler::PLVTypeIterator()) {
    CHECK_MESSAGE(plv_visited_map.find(plv_type) != plv_visited_map.end(),
                  "Iterator has PLV not in plv_vector.");
    plv_visited_map.at(plv_type) += 1;
  }
  // Check that each was visited only once.
  for (const auto [plv_type, visit_count] : plv_visited_map) {
    std::ignore = plv_type;
    CHECK_FALSE_MESSAGE(visit_count < 1, "One or more PLVs skipped by EnumIterator.");
    CHECK_FALSE_MESSAGE(visit_count > 1,
                        "One or more PLVs in visited more than once by EnumIterator.");
  }
}

#endif  // DOCTEST_LIBRARY_INCLUDED
