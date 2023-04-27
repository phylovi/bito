#pragma once

#include "../src/eigen_sugar.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED

void CheckVectorXdEquality(double value, const EigenVectorXd v, double tolerance) {
  for (Eigen::Index i = 0; i < v.size(); i++) {
    CHECK_LT(fabs(value - v[i]), tolerance);
  }
};

void CheckVectorXdEquality(const EigenVectorXd v1, const EigenVectorXd v2,
                           double tolerance) {
  CHECK_EQ(v1.size(), v2.size());
  for (Eigen::Index i = 0; i < v1.size(); i++) {
    double error = fabs(v1[i] - v2[i]);
    if (error > tolerance) {
      std::cerr << "CheckVectorXdEquality failed for index " << i << ": " << v1[i]
                << " vs " << v2[i] << std::endl;
    }
    CHECK_LT(error, tolerance);
  }
};

// Return the maximum absolute difference between any two entries in vector.
double VectorXdMaxError(const EigenVectorXd v1, const EigenVectorXd v2) {
  double max_error = 0.;
  Assert(v1.size() == v2.size(),
         "Cannot find max error of EigenVectorXd's of different sizes.");
  for (Eigen::Index i = 0; i < v1.size(); i++) {
    double error = fabs(v1[i] - v2[i]);
    if (error > max_error) {
      max_error = error;
    }
  }
  return max_error;
}

// Check if vectors are equal, within given tolerance for any two entries in vector.
bool VectorXdEquality(const EigenVectorXd v1, const EigenVectorXd v2,
                      double tolerance) {
  if (v1.size() != v2.size()) {
    return false;
  }
  for (Eigen::Index i = 0; i < v1.size(); i++) {
    double error = fabs(v1[i] - v2[i]);
    if (error > tolerance) {
      return false;
    }
  }
  return true;
};

void CheckVectorXdEqualityAfterSorting(const EigenVectorXdRef v1,
                                       const EigenVectorXdRef v2, double tolerance) {
  EigenVectorXd v1_sorted = v1;
  EigenVectorXd v2_sorted = v2;
  std::sort(v1_sorted.begin(), v1_sorted.end());
  std::sort(v2_sorted.begin(), v2_sorted.end());
  CheckVectorXdEquality(v1_sorted, v2_sorted, tolerance);
};

TEST_CASE(
    "Make sure that EigenVectorXdOfStdVectorDouble makes a new vector rather than "
    "wrapping data.") {
  std::vector<double> a = {1., 2., 3., 4.};
  EigenVectorXd b = EigenVectorXdOfStdVectorDouble(a);
  a[0] = 99;
  CHECK_EQ(b[0], 1.);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
