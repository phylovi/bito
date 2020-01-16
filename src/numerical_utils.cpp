// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "numerical_utils.hpp"
#include <iostream>

double NumericalUtils::LogSum(const EigenVectorXdRef vec) { return vec.redux(LogAdd); }

void NumericalUtils::ProbabilityNormalizeInLog(EigenVectorXdRef vec) {
  vec = vec.array() - LogSum(vec);
}

void NumericalUtils::Exponentiate(EigenVectorXdRef vec) { vec = vec.array().exp(); }
