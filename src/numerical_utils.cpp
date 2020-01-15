// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "numerical_utils.hpp"
#include <iostream>

using namespace std;

double NumericalUtils::LogSum(const EigenVectorXdRef vec) {
  double log_sum = DOUBLE_NEG_INF;
  for (auto val : vec) {
    log_sum = LogAdd(log_sum, val);
  }
  return log_sum;
}

void NumericalUtils::NormalizeInLog(EigenVectorXdRef vec) {
  double log_sum = LogSum(vec);
  for (auto &log_val : vec) {
    log_val -= log_sum;
  }
}

void NumericalUtils::Exponentiate(EigenVectorXdRef vec) {
  for (size_t i = 0; i < vec.size(); i++) {
    vec[i] = exp(vec[i]);
  }
}
