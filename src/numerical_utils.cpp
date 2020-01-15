//
//  numerical_utils.cpp
//  libsbn
//
//  Created by Seong-Hwan Jun on 1/13/20.
//

#include "numerical_utils.hpp"
#include <iostream>

using namespace std;

double NumericalUtils::LogAdd(double x, double y) {
  // make x the max
  if (y > x) {
    double temp = x;
    x = y;
    y = temp;
  }
  // now x is bigger
  if (x == DOUBLE_NEG_INF) {
    return x;
  }
  double negDiff = y - x;
  if (negDiff < -20) {
    return x;
  }
  return x + log(1.0 + exp(negDiff));
}

double NumericalUtils::LogSum(const EigenVectorXdRef vec) {
  double log_sum = DOUBLE_NEG_INF;
  for (auto val : vec) {
    log_sum = LogAdd(log_sum, val);
  }
  return log_sum;
}

void NumericalUtils::LogNormalize(EigenVectorXdRef vec) {
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
