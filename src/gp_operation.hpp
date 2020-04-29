// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef GP_OPERATION_HPP_
#define GP_OPERATION_HPP_

#include <iostream>
#include <variant>
#include <vector>
#include "sugar.hpp"

using StringSizePairVector = std::vector<std::pair<std::string, size_t>>;

namespace GPOperations {

struct Zero {
  size_t idx;

  StringSizePairVector guts() const { return {{"idx", idx}}; }
};

struct WeightedSumAccumulate {
  size_t dest_idx;
  size_t weight_idx;
  size_t src_idx;

  StringSizePairVector guts() const {
    return {{"dest_idx", dest_idx}, {"weight_idx", weight_idx}, {"src_idx", src_idx}};
  }
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

using GPOperationVector = std::vector<GPOperation>;

struct GPOperationOstream {
  std::ostream& os_;

  GPOperationOstream(std::ostream& os) : os_{os} {}

  void operator()(const GPOperations::Zero& operation) {
    os_ << "Zero" << operation.guts();
  }
  void operator()(const GPOperations::WeightedSumAccumulate& operation) {
    os_ << "WeightedSumAccumulate" << operation.guts();
  }
};

std::ostream& operator<<(std::ostream& os, GPOperation const& v) {
  std::visit(GPOperationOstream{os}, v);
  return os;
}


#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("GPOperation") {
  // TODO
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // GP_OPERATION_HPP_
