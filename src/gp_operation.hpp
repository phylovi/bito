// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef GP_OPERATION_HPP_
#define GP_OPERATION_HPP_

#include <iostream>
#include <variant>
#include <vector>
#include "sugar.hpp"

namespace GPOperations {

struct Zero {
  size_t idx;
};

struct WeightedSumAccumulate {
  size_t dest_idx;
  size_t weight_idx;
  size_t src_idx;
};

struct Multiply {
  size_t dest_idx;
  size_t src1_idx;
  size_t src2_idx;
};

struct Likelihood {
  size_t dest_idx;
  size_t src1_idx;
  size_t src2_idx;
};

struct EvolveRootward {
  size_t dest_idx;
  size_t src_idx;
  size_t branch_length_idx;
};

struct EvolveLeafward {
  size_t dest_idx;
  size_t src_idx;
  size_t branch_length_idx;
};

struct OptimizeRootward {
  size_t plv_dest_idx;
  size_t leafward_idx;
  size_t rootward_idx;
  size_t branch_length_idx;
};

struct OptimizeLeafward {
  size_t plv_dest_idx;
  size_t leafward_idx;
  size_t rootward_idx;
  size_t branch_length_idx;
};

struct UpdateSBNProbabilities {
  size_t start_idx;
  size_t stop_idx;
};
}  // namespace GPOperations

using GPOperation =
    std::variant<GPOperations::Zero, GPOperations::WeightedSumAccumulate>;
//                 GPOperations::Multiply, GPOperations::Likelihood,
//                 GPOperations::EvolveRootward, GPOperations::EvolveLeafward,
//                 GPOperations::OptimizeRootward, GPOperations::OptimizeLeafward,
//                 GPOperations::UpdateSBNProbabilities>;

struct GPToMap {
  StringSizeMap operator()(const GPOperations::Zero& operation) {
    return {{"idx", operation.idx}};
  }
  StringSizeMap operator()(const GPOperations::WeightedSumAccumulate& operation) {
    return {{"dest_idx", operation.dest_idx},
            {"weight_idx", operation.weight_idx},
            {"src_idx", operation.src_idx}};
  }
};

struct GPOperationOstream {
  std::ostream& os_;

  GPOperationOstream(std::ostream& os) : os_{os} {}

  void operator()(const GPOperations::Zero& operation) {
    os_ << "Zero" << std::visit(GPToMap{}, operation);
  }
  void operator()(const GPOperations::WeightedSumAccumulate& operation) {
    os_ << "WeightedSumAccumulate" << std::visit(GPToMap{}, operation);
  }
};

std::ostream& operator<<(std::ostream& os, GPOperation const& v) {
  std::visit(GPOperationOstream{os}, v);
  return os;
}

// struct GPOperationOstream {
//   std::ostream& os_;
//
//   GPOperationOstream(std::ostream& os) : os_{os} {}
//
//   template <typename... Args>
//   void FoldPushBack(std::vector<size_t>& v, Args&&... args) {
//     (v.push_back(args), ...);
//   }
//
//   template <typename... Args>
//   std::vector<size_t> MakeVector(Args&&... args) {
//     std::vector<size_t> v;
//     FoldPushBack(v, &args...);
//     return v;
//   }
//
//   template <typename... Args>
//   void FoldPrint(Args&&... args) {
//     (os_ << ... << std::forward<Args>(args)) << '\n';
//   }
//
//   void operator()(const GPOperations::Zero& operation) {
//     FoldPrint("Zero", operation.idx);
//   }
//   void operator()(const GPOperations::WeightedSumAccumulate& operation) {
//     auto v = MakeVector(operation.dest_idx, operation.weight_idx, operation.src_idx);
//     // std::vector<size_t> v{operation.dest_idx, operation.weight_idx,
//     // operation.src_idx};
//     os_ << "WeightedSumAccumulate" << v;
//     // FoldPrint("WeightedSumAccumulate", operation.dest_idx, operation.weight_idx,
//     //           operation.src_idx);
//   }
// };
//
// std::ostream& operator<<(std::ostream& os, GPOperation const& v) {
//   std::visit(GPOperationOstream{os}, v);
//   return os;
// }

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("TODO") {}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // GP_OPERATION_HPP_
