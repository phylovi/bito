// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <fenv.h>

#include <limits>

#include "eigen_sugar.hpp"
#include "sugar.hpp"

constexpr double DOUBLE_INF = std::numeric_limits<double>::infinity();
constexpr double DOUBLE_NEG_INF = -std::numeric_limits<double>::infinity();
constexpr double EPS = std::numeric_limits<double>::epsilon();
constexpr size_t OUT_OF_SAMPLE_IDX = std::numeric_limits<size_t>::max();
constexpr double DOUBLE_NAN = std::numeric_limits<double>::quiet_NaN();
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
inline double DOUBLE_MINIMUM = std::numeric_limits<double>::lowest() * ERR_TOLERANCE;

constexpr auto FE_OVER_AND_UNDER_FLOW_EXCEPT = FE_OVERFLOW | FE_UNDERFLOW;

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
