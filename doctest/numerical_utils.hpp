#pragma once

#include "../src/numerical_utils.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("NumericalUtils") {
  double log_x = log(2);
  double log_y = log(3);
  double log_sum = NumericalUtils::LogAdd(log_x, log_y);
  CHECK_LT(fabs(log_sum - 1.609438), 1e-5);

  EigenVectorXd log_vec(10);
  double log_sum2 = DOUBLE_NEG_INF;
  for (Eigen::Index i = 0; i < log_vec.size(); i++) {
    log_vec(i) = log(i + 1);
    log_sum2 = NumericalUtils::LogAdd(log_sum2, log_vec(i));
  }
  log_sum = NumericalUtils::LogSum(log_vec);
  CHECK_LT(fabs(log_sum - 4.007333), 1e-5);
  CHECK_LT(fabs(log_sum2 - 4.007333), 1e-5);

  NumericalUtils::ProbabilityNormalizeInLog(log_vec);
  for (Eigen::Index i = 0; i < log_vec.size(); i++) {
    CHECK_LT(fabs(log_vec(i) - (log(i + 1) - log_sum)), 1e-5);
  }

  NumericalUtils::Exponentiate(log_vec);
  double sum = 0.0;
  for (Eigen::Index i = 0; i < log_vec.size(); i++) {
    sum += log_vec(i);
  }
  CHECK_LT(fabs(sum - 1), 1e-5);

  // Here we use volatile to avoid GCC optimizing away the variable.
  volatile double d = 4.;
  std::ignore = d;
  d /= 0.;
  auto fp_description = NumericalUtils::DescribeFloatingPointEnvironmentExceptions();
  CHECK_EQ(*fp_description,
           "The following floating point problems have been encountered: FE_DIVBYZERO");
}
#endif  // DOCTEST_LIBRARY_INCLUDED
