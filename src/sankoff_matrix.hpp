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

#ifdef DOCTEST_LIBRARY_INCLUDED
enum Nucleotides { A, C, G, T };

TEST_CASE("SankoffMatrix: Testing SankoffMatrix Getter/Setter Methods") {
  auto cm = SankoffMatrix();
  CHECK_LT(fabs(cm.GetCost(0, 0) - 0.), 1e-10);
  CHECK_LT(fabs(cm.GetCost(G, C) - 1.), 1e-10);

  cm.UpdateMatrix(A, G, 3.);
  auto test_matrix = Eigen::Matrix<double, 4, 4>();
  test_matrix << 0., 1., 3., 1., 1., 0., 1., 1., 1., 1., 0., 1., 1., 1., 1., 0.;
  CHECK(test_matrix.isApprox(cm.GetMatrix()));
  CHECK_LT(fabs(cm.GetCost(A, G) - 3.), 1e-10);
  CHECK_THROWS(cm.UpdateMatrix(G, G, 3.));
}

TEST_CASE("SankoffMatrix: Create SankoffMatrix from given cost matrix") {
  auto costs = Eigen::Matrix<double, 4, 4>();
  costs << 0., 2.5, 1., 2.5, 2.5, 0., 2.5, 1., 1., 2.5, 0., 2.5, 2.5, 1., 2.5, 0.;
  auto cm = SankoffMatrix(costs);
  CHECK_LT(fabs(cm.GetCost(A, C) - 2.5), 1e-10);

  auto costs_invalid = Eigen::Matrix<double, 4, 4>();
  // non-zero values on diagonal
  costs_invalid << 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.,
      16.;
  CHECK_THROWS(new SankoffMatrix(costs_invalid));
}

#endif  // DOCTEST_LIBRARY_INCLUDED
