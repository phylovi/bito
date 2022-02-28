// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// PLVHandler is used by GPOperations to get index location of PLVs stored on the
// GPEngine.

#pragma once

class PLVHandler {
 public:
  // We store 6 PLVs per subsplit, and index them according to this enum. (See
  // GPEngine::plvs_).
  static inline const size_t plv_count = 6;
  enum class PLVType : size_t {
    P,          // p(s)
    PHatRight,  // phat(s_right)
    PHatLeft,   // phat(s_left)
    RHat,       // rhat(s_right) = rhat(s_left)
    RRight,     // r(s_right)
    RLeft       // r(s_left)
  };

  // Get the `GPEngine::plvs_` index of given node's given PLV type from DAG with
  // node_count.
  static size_t GetPLVIndex(const PLVType plv_type, const size_t node_count,
                            const size_t node_idx) {
    size_t plv_type_idx;
    switch (plv_type) {
      case PLVType::P:
        plv_type_idx = 0;
        break;
      case PLVType::PHatRight:
        plv_type_idx = 1;
        break;
      case PLVType::PHatLeft:
        plv_type_idx = 2;
        break;
      case PLVType::RHat:
        plv_type_idx = 3;
        break;
      case PLVType::RRight:
        plv_type_idx = 4;
        break;
      case PLVType::RLeft:
        plv_type_idx = 5;
        break;
      default:
        Failwith("Invalid PLV index requested.");
    }
    return (plv_type_idx * node_count) + node_idx;
  };

  static size_t GetPLVIndex(const PLVType plv_type, const SubsplitDAG &dag,
                            const size_t node_idx) {
    return GetPLVIndex(plv_type, dag.NodeCountWithoutDAGRoot(), node_idx);
  };

  static PLVType RPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::RLeft : PLVType::RRight;
  };
};
