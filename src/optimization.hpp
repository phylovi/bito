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
std::tuple<T, T, std::tuple<std::vector<T>, std::vector<T>, std::vector<T>>>
BrentMinimize(F f, T guess, T min, T max, int significant_digits, size_t max_iter) {
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
  T f_prime_x;

  // golden ratio, don't need too much precision here!
  static const T golden = 0.3819660f;

  std::vector<T> path_x;
  std::vector<T> path_fx;
  std::vector<T> path_fprimex;


  w = v = x = guess;
  fw = fv = fx = f(x).first;
  delta2 = delta = 0;

  path_x.push_back(x);
  path_fx.push_back(fx);
  path_fprimex.push_back(f(x).second);

  size_t count = max_iter;

  do {
    // Check current value
    f_prime_x = f(x).second;
    if (fabs(f_prime_x) < 1e-5) {
      break;
    }

    // get midpoint
    mid = (min + max) / 2;

    // ** Convergence Test:
    // work out if we're done already:
    fract1 = tolerance * fabs(x) + tolerance / 4;
    fract2 = 2 * fract1;
    if (fabs(x - mid) <= (fract2 - (max - min) / 2)) {
      break;
    }

    // ** Parabolic Fit (variant of Inverse Quadratic Interpolation?):
    bool use_bisection = true;  // only triggers if IQI fails.
    if (fabs(delta2) > fract1) {
      // Try and construct a parabolic fit:
      T r = (x - w) * (fx - fv);
      T q = (x - v) * (fx - fw);
      T p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) p = -p;
      q = fabs(q);
      T td = delta2;
      delta2 = delta;
      // Determine whether a parabolic step is acceptable or not:
      // must fail all three threshold tests to be accepted.
      if (((fabs(p) >= fabs(q * td / 2)) == false) && ((p <= q * (min - x)) == false) &&
          ((p >= q * (max - x)) == false)) {
        // whew, parabolic fit:
        delta = p / q;
        u = x + delta;
        if (((u - min) < fract2) || ((max - u) < fract2)) {
          delta = (mid - x) < 0 ? (T)-fabs(fract1) : (T)fabs(fract1);
        }
        // parabolic fit was a success, so don't need bisection.
        use_bisection = false;
      }
    }

    // ** Golden Bisection Method (this is an optimization of traditional Bisection
    // Method)
    if (use_bisection) {
      // golden section:
      delta2 = (x >= mid) ? min - x : max - x;
      delta = golden * delta2;
    }

    // ** Update current position:
    u = (fabs(delta) >= fract1)
            ? T(x + delta)
            : (delta > 0 ? T(x + fabs(fract1)) : T(x - fabs(fract1)));
    fu = f(u).first;
    if (fu <= fx) {
      // good new point is an improvement!
      // update brackets (previous guess becomes the new outer bracket):
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

    path_x.push_back(x);
    path_fx.push_back(fx);
    path_fprimex.push_back(f(x).second);

  } while (--count);  // countdown until max iterations.

  max_iter -= count;

  return std::make_tuple(x, fx, std::make_tuple(path_x, path_fx, path_fprimex));
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

//
DoublePair LogSpaceGradientAscent(std::function<DoublePair(double)> f_and_f_prime,
                                  double x, const double tolerance,
                                  const double log_space_step_size, const double min_x,
                                  const size_t max_iter) {
  size_t iter_idx = 0;
  while (true) {
    double y = log(x);
    auto [f_x, f_prime_x] = f_and_f_prime(x);
    double log_space_grad = x * f_prime_x;
    const double new_y = y + log_space_grad * log_space_step_size;
    const double new_x = exp(new_y);
    x = std::max(new_x, min_x);
    if (fabs(f_prime_x) < fabs(f_x) * tolerance || iter_idx >= max_iter) {
      return {x, f_x};
    }
    ++iter_idx;
  }
}

// Modifying the output so that we can plot the optimization path
// TODO: Change return back to DoublePair before merging to main.
std::tuple<double, double,
           std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>
NewtonRaphsonOptimization(
    std::function<std::tuple<double, double, double>(double)> f_and_derivatives,
    double x, const double tolerance, const double epsilon, const double min_x,
    const double max_x, const size_t max_iter) {
  size_t iter_idx = 0;
  double new_x, delta;
  double damp_const = 0.005;

  std::vector<double> path_x;
  std::vector<double> path_fx;
  std::vector<double> path_fprimex;

  while (true) {
    auto [f_x, f_prime_x, f_double_prime_x] = f_and_derivatives(x);
    path_x.push_back(x);
    path_fx.push_back(f_x);
    path_fprimex.push_back(f_prime_x);

    if (fabs(f_double_prime_x) < epsilon) {
      return std::make_tuple(x, f_x, std::make_tuple(path_x, path_fx, path_fprimex));
    }
    // This method below produces reasonable estimates on single tree DS data
    //
    new_x = x - damp_const * f_prime_x / f_double_prime_x;

    if (new_x <= min_x) {
      new_x = x - 0.5 * (x - min_x);
      damp_const *= 10.;
    } else if (new_x >= max_x) {
      new_x = x - 0.5 * (x - max_x);
      damp_const *= 10.;
    } else {
      // damp_const *= 0.5;
    }

    delta = fabs(x - new_x);

    if (delta < tolerance || iter_idx == max_iter) {
      return std::make_tuple(x, f_x, std::make_tuple(path_x, path_fx, path_fprimex));
    }

    x = new_x;
    ++iter_idx;
  }
}

DoublePair LogSpaceNewtonRaphsonOptimization(
    std::function<std::tuple<double, double, double>(double)> f_and_derivatives,
    double x, const double tolerance, const double epsilon, const double min_x,
    const size_t max_iter) {
  size_t iter_idx = 0;
  while (true) {
    auto [f_x, f_prime_x, f_double_prime_x] = f_and_derivatives(x);

    double y = log(x);
    double f_prime_y = x * f_prime_x;
    double f_double_prime_y = f_prime_y + pow(x, 2) * f_double_prime_x;

    const double new_y = y - f_prime_y / f_double_prime_y;
    double new_x = fmax(exp(new_y), min_x);
    double delta = fabs(x - new_x);

    if (delta < tolerance || fabs(f_double_prime_y) < epsilon || iter_idx == max_iter) {
      return {x, f_x};
    } else {
      x = new_x;
    }
    ++iter_idx;
  }
}

}  // namespace Optimization
