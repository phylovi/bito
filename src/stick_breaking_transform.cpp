// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This code closely follows
// https://mc-stan.org/docs/2_26/reference-manual/simplex-transform-section.html
// and so we follow their notation.

#include "stick_breaking_transform.hpp"

inline double inverse_logit(const double y) { return 1.0 / (1 + std::exp(-y)); }

inline double logit(const double x) { return std::log(x / (1.0 - x)); }

inline double log1p_exp(double a) {
  // prevents underflow
  if (a > 0.0) return a + std::log1p(exp(-a));
  return std::log1p(exp(a));
}

EigenVectorXd StickBreakingTransform::operator()(EigenVectorXd const& y) const {
  size_t K = y.size() + 1;
  EigenVectorXd x(K);
  double stick = 1.0;
  for (size_t k = 0; k < K - 1; k++) {
    double z = inverse_logit(y[k] - std::log(K - k - 1));
    x[k] = stick * z;
    stick -= x[k];
  }
  x[K - 1] = stick;
  return x;
}

EigenVectorXd StickBreakingTransform::inverse(const EigenVectorXd& x) const {
  size_t K = x.size();
  double sum = 0;
  EigenVectorXd y(K - 1);
  for (size_t k = 0; k < K - 1; k++) {
    double z = x[k] / (1.0 - sum);
    y[k] = logit(z) + std::log(K - k - 1);
    sum += x[k];
  }
  return y;
}

double StickBreakingTransform::log_abs_det_jacobian(const EigenVectorXd& x,
                                                    const EigenVectorXd& y) const {
  size_t K = x.size();
  double log_prob = 0.0;
  double stick = 1.0;
  for (size_t k = 0; k < K - 1; k++) {
    double adj_y_k = y[k] - log(K - k - 1);
    log_prob += std::log(stick);
    log_prob -= log1p_exp(-adj_y_k);
    log_prob -= log1p_exp(adj_y_k);
    stick -= x[k];
  }
  return log_prob;
}
