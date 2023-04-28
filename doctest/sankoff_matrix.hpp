#pragma once

#include "../src/sankoff_matrix.hpp"

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
