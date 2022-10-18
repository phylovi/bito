// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "optimization.hpp"

double Optimization::GradientAscent(std::function<DoublePair(double)> f_and_f_prime,
                                    double x, const int significant_digits,
                                    const double step_size, const double min_x,
                                    const size_t max_iter) {
  double tolerance = std::pow(10, -significant_digits);
  size_t iter_idx = 0;
  while (true) {
    auto [f_x, f_prime_x] = f_and_f_prime(x);
    const double new_x = x + f_prime_x * step_size;
    x = std::max(new_x, min_x);
    if (fabs(f_prime_x) < fabs(f_x) * tolerance || iter_idx >= max_iter) {
      return x;
    }
    ++iter_idx;
  }
}

double Optimization::LogSpaceGradientAscent(
    std::function<DoublePair(double)> f_and_f_prime, double x,
    const int significant_digits, const double log_space_step_size, const double min_x,
    const size_t max_iter) {
  double tolerance = static_cast<double>(std::pow(10, -significant_digits));
  size_t iter_idx = 0;
  while (true) {
    double y = log(x);
    auto [f_x, f_prime_x] = f_and_f_prime(x);
    double log_space_grad = x * f_prime_x;
    const double new_y = y + log_space_grad * log_space_step_size;
    const double new_x = exp(new_y);
    x = std::max(new_x, min_x);
    if (fabs(f_prime_x) < fabs(f_x) * tolerance || iter_idx >= max_iter) {
      return x;
    }
    ++iter_idx;
  }
}

double Optimization::NewtonRaphsonOptimization(
    std::function<std::tuple<double, double, double>(double)> f_and_derivatives,
    double x, const int significant_digits, const double epsilon, const double min_x,
    const double max_x, const size_t max_iter) {
  double tolerance = pow(10, -significant_digits);
  size_t iter_idx = 0;
  double new_x, delta;
  double min = min_x;
  double max = max_x;

  while (true) {
    auto [f_x, f_prime_x, f_double_prime_x] = f_and_derivatives(x);

    if (fabs(f_double_prime_x) < epsilon) {
      return x;
    }
    new_x = x - f_prime_x / f_double_prime_x;

    if (new_x < min_x) {
      new_x = x - 0.5 * (x - min);
    }
    if (new_x > max_x) {
      new_x = x - 0.5 * (x - max);
    }

    delta = fabs(x - new_x);

    if (delta < tolerance || fabs(f_prime_x) < fabs(f_x) * tolerance ||
        iter_idx == max_iter) {
      return x;
    }

    x = new_x;
    ++iter_idx;
  }
}
