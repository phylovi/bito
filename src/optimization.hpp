// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <cmath>
#include <functional>
#include <numeric>
#include <utility>

#include "sugar.hpp"

namespace Optimization {

// Copied from https://www.boost.org/doc/libs/1_73_0/boost/math/tools/minima.hpp
template <class F, class T>
std::pair<T, T> BrentMinimize(F f, T min, T max, int significant_digits,
                              size_t max_iter) {
  T tolerance = static_cast<T>(ldexp(1.0, 1 - significant_digits));
  T x;               // minima so far
  T w;               // second best point
  T v;               // previous value of w
  T u;               // most recent evaluation point
  T delta;           // The distance moved in the last step
  T delta2;          // The distance moved in the step before last
  T fu, fv, fw, fx;  // function evaluations at u, v, w, x
  T mid;             // midpoint of min and max
  T fract1, fract2;  // minimal relative movement in x

  static const T golden =
      0.3819660f;  // golden ratio, don't need too much precision here!

  x = w = v = max;
  fw = fv = fx = f(x);
  delta2 = delta = 0;

  size_t count = max_iter;

  do {
    // get midpoint
    mid = (min + max) / 2;
    // work out if we're done already:
    fract1 = tolerance * fabs(x) + tolerance / 4;
    fract2 = 2 * fract1;
    if (fabs(x - mid) <= (fract2 - (max - min) / 2)) {
      break;
    }

    if (fabs(delta2) > fract1) {
      // try and construct a parabolic fit:
      T r = (x - w) * (fx - fv);
      T q = (x - v) * (fx - fw);
      T p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) p = -p;
      q = fabs(q);
      T td = delta2;
      delta2 = delta;
      // determine whether a parabolic step is acceptable or not:
      if ((fabs(p) >= fabs(q * td / 2)) || (p <= q * (min - x)) ||
          (p >= q * (max - x))) {
        // nope, try golden section instead
        delta2 = (x >= mid) ? min - x : max - x;
        delta = golden * delta2;
      } else {
        // whew, parabolic fit:
        delta = p / q;
        u = x + delta;
        if (((u - min) < fract2) || ((max - u) < fract2))
          delta = static_cast<T>((mid - x) < 0 ? -fabs(fract1) : fabs(fract1));
      }
    } else {
      // golden section:
      delta2 = (x >= mid) ? min - x : max - x;
      delta = golden * delta2;
    }
    // update current position:
    u = (fabs(delta) >= fract1)
            ? T(x + delta)
            : (delta > 0 ? T(x + fabs(fract1)) : T(x - fabs(fract1)));
    fu = f(u);
    if (fu <= fx) {
      // good new point is an improvement!
      // update brackets:
      if (u >= x)
        min = x;
      else
        max = x;
      // update control points:
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    } else {
      // Oh dear, point u is worse than what we have already,
      // even so it *must* be better than one of our endpoints:
      if (u < x)
        min = u;
      else
        max = u;
      if ((fu <= fw) || (w == x)) {
        // however it is at least second best:
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } else if ((fu <= fv) || (v == x) || (v == w)) {
        // third best:
        v = u;
        fv = fu;
      }
    }

  } while (--count);

  max_iter -= count;

  return std::make_pair(x, fx);
}

DoublePair GradientAscent(std::function<DoublePair(double)> f_and_f_prime, double x,
                          const double tolerance, const double step_size,
                          const double min_x, const size_t max_iter) {
  size_t iter_idx = 0;
  while (true) {
    auto [f_x, f_prime_x] = f_and_f_prime(x);
    const double new_x = x + f_prime_x * step_size;
    x = std::max(new_x, min_x);
    if (fabs(f_prime_x) < fabs(f_x) * tolerance || iter_idx >= max_iter) {
      return {x, f_x};
    }
    ++iter_idx;
  }
}

DoublePair LogSpaceGradientAscent(std::function<DoublePair(double)> f_and_f_prime,
                                  double x, const double tolerance,
                                  const double step_size, const double min_x,
                                  const size_t max_iter) {
  size_t iter_idx = 0;
  while (true) {
    double log_x = log(x);
    auto [f_x, f_prime_x] = f_and_f_prime(x);
    auto [f_log_x, f_prime_log_x] = f_and_f_prime(log_x);
    double log_space_grad = (1. / x) * f_prime_log_x;
    const double new_log_x = log_x + log_space_grad * step_size;
    x = std::max(exp(new_log_x), min_x);
    if (fabs(log_space_grad) < fabs(f_log_x) * tolerance || iter_idx >= max_iter) {
      return {x, f_x};
    }
    ++iter_idx;
  }
}

DoublePair AdaptiveGradientAscent(
    std::function<DoublePair(double)> f_and_f_prime, double x, const double tolerance,
    const double default_step_size, const double uniformity_bound,
    const size_t monotonicity_const, const double threshold_const,
    const double rescaling_constant, const size_t max_iter) {
  size_t iter_idx = 0;
  double alpha = 1 / default_step_size;
  std::vector<double> f_log_values(monotonicity_const, 0.);

  while (true) {
    double log_x = log(x);
    auto [f_x, f_prime_x] = f_and_f_prime(x);
    auto [f_log_x, f_prime_log_x] = f_and_f_prime(log_x);
    double log_space_grad = (1. / x) * f_prime_log_x;
    double log_space_grad_sq = std::pow(log_space_grad, 2);

    int f_log_value_index = static_cast<int>(iter_idx % monotonicity_const);
    f_log_values.at(f_log_value_index) = f_log_x;
    // int last_element_index =
    //    std::min(static_cast<int>(iter_idx + 1),
    //    static_cast<int>(monotonicity_const));
    double min_f_log = *std::min_element(f_log_values.begin(), f_log_values.end());

    if (alpha <= uniformity_bound || alpha >= (1 / uniformity_bound)) {
      alpha = 1 / default_step_size;
    }
    double lambda_step_size = 1 / alpha;

    double candidate_log_x = std::max(log_x + log_space_grad * lambda_step_size, log_x);
    auto [f_candidate_log_x, f_prime_candidate_log_x] = f_and_f_prime(candidate_log_x);
    double acceptance_threshold =
        threshold_const * lambda_step_size * log_space_grad_sq;

    while (f_candidate_log_x < min_f_log + acceptance_threshold) {
      lambda_step_size = lambda_step_size * rescaling_constant;
      candidate_log_x = std::max(log_x + log_space_grad * lambda_step_size, log_x);
      auto [f_candidate_log_x, f_prime_candidate_log_x] =
          f_and_f_prime(candidate_log_x);
      acceptance_threshold =
          threshold_const * lambda_step_size * log_space_grad_sq;
    }

    x = exp(candidate_log_x);
    alpha = (-log_space_grad *
             (f_prime_candidate_log_x * exp(-candidate_log_x) - log_space_grad)) /
            (lambda_step_size * log_space_grad_sq);

    if (fabs(log_space_grad) < fabs(f_log_x) * tolerance || iter_idx >= max_iter) {
      return {x, f_x};
    }
    // TODO: should overwrite f_x and f_prime_x with the candidate values
    // so that they are not unnecessarily recomputed at the top of the loop
    ++iter_idx;
  }
}
}  // namespace Optimization
