// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_NUMERICAL_UTILS_HPP_
#define SRC_NUMERICAL_UTILS_HPP_

#include <fenv.h>
#include <limits>
#include <string>

#include "eigen_sugar.hpp"
#include "sugar.hpp"

constexpr double DOUBLE_INF = std::numeric_limits<double>::infinity();
constexpr double DOUBLE_NEG_INF = -std::numeric_limits<double>::infinity();
constexpr double EPS = std::numeric_limits<double>::epsilon();
// It turns out that log isn't constexpr for silly reasons, so we use inline instead.
inline double LOG_EPS = log(EPS);
constexpr double ERR_TOLERANCE = 1e-10;
// DOUBLE_MINIMUM defines de facto minimum value for double to deal with
// potential overflow resulting from summing of large number of log
// probabilities.
// Note: using std::numeric_limits<double>::lowest() may result in numerical
// instability, especially if other operations are to be performed using it.
// This is why we are using a value that is slightly larger to denote
// the lowest double value that we will consider.
inline double DOUBLE_MINIMUM = std::numeric_limits<double>::lowest()*ERR_TOLERANCE;

namespace NumericalUtils {

// Return log(exp(x) + exp(y)).
inline double LogAdd(double x, double y) {
  // See:
  // https://github.com/alexandrebouchard/bayonet/blob/master/src/main/java/bayonet/math/NumericalUtils.java#L59
  // Make x the max.
  if (y > x) {
    double temp = x;
    x = y;
    y = temp;
  }
  if (x == DOUBLE_NEG_INF) {
    return x;
  }
  double neg_diff = y - x;
  if (neg_diff < LOG_EPS) {
    return x;
  }
  return x + log(1.0 + exp(neg_diff));
}

// Return log(sum_i exp(vec(i))).
double LogSum(const EigenVectorXdRef vec);
// Returns a vector with the i-th entry given by LogAdd(vec1(i), vec2(i))
EigenVectorXd LogAddVectors(const EigenVectorXdRef vec1, const EigenVectorXdRef vec2);
// Normalize the entries of vec such that they become logs of probabilities:
// vec(i) = vec(i) - LogSum(vec).
void ProbabilityNormalizeInLog(EigenVectorXdRef vec);
// Exponentiate vec in place: vec(i) = exp(vec(i))
void Exponentiate(EigenVectorXdRef vec);

// This code concerns the floating-point environment (FE).
// Note that a FE "exception" is not a C++ exception, it's just a signal that a
// problem has happened. See
// https://en.cppreference.com/w/c/numeric/fenv/FE_exceptions

// If there is a worrying FE exception, get a string describing it.
std::optional<std::string> DescribeFloatingPointEnvironmentExceptions();
// If there is a worrying FE exception, report it and clear the record of there being
// any exception.
void ReportFloatingPointEnvironmentExceptions(std::string context = "");

}  // namespace NumericalUtils

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("NumericalUtils") {
  double log_x = log(2);
  double log_y = log(3);
  double log_sum = NumericalUtils::LogAdd(log_x, log_y);
  CHECK_LT(fabs(log_sum - 1.609438), 1e-5);

  EigenVectorXd log_vec(10);
  double log_sum2 = DOUBLE_NEG_INF;
  for (size_t i = 0; i < log_vec.size(); i++) {
    log_vec(i) = log(i + 1);
    log_sum2 = NumericalUtils::LogAdd(log_sum2, log_vec(i));
  }
  log_sum = NumericalUtils::LogSum(log_vec);
  CHECK_LT(fabs(log_sum - 4.007333), 1e-5);
  CHECK_LT(fabs(log_sum2 - 4.007333), 1e-5);

  NumericalUtils::ProbabilityNormalizeInLog(log_vec);
  for (size_t i = 0; i < log_vec.size(); i++) {
    CHECK_LT(fabs(log_vec(i) - (log(i + 1) - log_sum)), 1e-5);
  }

  NumericalUtils::Exponentiate(log_vec);
  double sum = 0.0;
  for (size_t i = 0; i < log_vec.size(); i++) {
    sum += log_vec(i);
  }
  CHECK_LT(fabs(sum - 1), 1e-5);

  // Here we use volatile to avoid GCC optimizing away the variable.
  volatile double d = 4.;
  d /= 0.;
  auto fp_description = NumericalUtils::DescribeFloatingPointEnvironmentExceptions();
  CHECK_EQ(*fp_description,
           "The following floating point problems have been encountered: FE_DIVBYZERO");
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif /* SRC_NUMERICAL_UTILS_HPP_ */
