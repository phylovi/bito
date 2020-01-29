// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "numerical_utils.hpp"
#include <iostream>

double NumericalUtils::LogSum(const EigenVectorXdRef vec) { return vec.redux(LogAdd); }

EigenVectorXd NumericalUtils::LogAddVectors(const EigenVectorXdRef vec1,
                                            const EigenVectorXdRef vec2) {
  EigenVectorXd result(vec1.size());
  Assert(vec1.size() == vec2.size(), "Two vectors must have the same length.");
  std::transform(vec1.begin(), vec1.end(), vec2.begin(), result.begin(), LogAdd);
  return result;
}

void NumericalUtils::ProbabilityNormalizeInLog(EigenVectorXdRef vec) {
  vec = vec.array() - LogSum(vec);
}

void NumericalUtils::Exponentiate(EigenVectorXdRef vec) { vec = vec.array().exp(); }

// This is any FE exception except for FE_INEXACT, which happens all the time.
constexpr auto FE_WORRYING_EXCEPT =
    FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW;

std::optional<std::string>
NumericalUtils::DescribeFloatingPointEnvironmentExceptions() {
  if (fetestexcept(FE_WORRYING_EXCEPT)) {
    std::string warning_string(
        "The following floating point problems have been encountered:");
    if (fetestexcept(FE_DIVBYZERO)) {
      warning_string.append(" FE_DIVBYZERO");
    }
    if (fetestexcept(FE_INVALID)) {
      warning_string.append(" FE_INVALID");
    }
    if (fetestexcept(FE_OVERFLOW)) {
      warning_string.append(" FE_OVERFLOW");
    }
    if (fetestexcept(FE_UNDERFLOW)) {
      warning_string.append(" FE_UNDERFLOW");
    }
    return warning_string;
  }  // else
  return std::nullopt;
}

void NumericalUtils::ReportFloatingPointEnvironmentExceptions(std::string context) {
  auto warning_string = DescribeFloatingPointEnvironmentExceptions();
  if (warning_string) {
    std::cout << context << " " << *warning_string << std::endl;
    feclearexcept(FE_ALL_EXCEPT);
  }
}
