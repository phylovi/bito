// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Cost Matrix for Sankoff Algorithm
// cost_matrix[i][j] is the cost of mutating from parent state i to child state j

#pragma once

#include "eigen_sugar.hpp"
#include "sugar.hpp"

using CostMatrix = Eigen::Matrix<double, 4, 4>;

class SankoffMatrix {
 public:
  static const size_t state_count = 4;

  // default matrix is all ones except for zeroes on the diagonal
  SankoffMatrix() : cost_matrix_(CostMatrix()) {
    cost_matrix_.setOnes();
    for (size_t state = 0; state < state_count; state++) {
      cost_matrix_(state, state) = 0.;
    }
  };
  SankoffMatrix(CostMatrix cost_matrix) {
    for (size_t state = 0; state < state_count; state++) {
      Assert(fabs(cost_matrix(state, state) - 0.) < 1e-10,
             "Diagnonal of cost matrix should be 0.");
    }
    cost_matrix_ = std::move(cost_matrix);
  }

  void UpdateMatrix(size_t parent_state, size_t child_state, double cost) {
    Assert(parent_state != child_state,
           "Mutating from state " + std::to_string(parent_state) + " to state " +
               std::to_string(child_state) + " should have cost of 0.");
    cost_matrix_(parent_state, child_state) = cost;
  };
  double GetCost(size_t parent_state, size_t child_state) {
    return cost_matrix_(parent_state, child_state);
  };
  CostMatrix GetMatrix() { return cost_matrix_; };

 private:
  CostMatrix cost_matrix_;
};
