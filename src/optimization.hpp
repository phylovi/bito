// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <cmath>
#include <functional>
#include <numeric>
#include <utility>

#include "sugar.hpp"

// We tested both gradient and gradient-free based methods for branch length
// optimization using generalized pruning. We found that gradient-based, particularly
// Newton's method, works well on small test cases, but fails compared to gradient-free
// methods in real data. We tried a "hybrid" approach, where ran Brent for the first
// round or 2 of optimization (on the full DAG), and then handed off to Newton's method,
// but these didn't show improvements over Newton's method.
//
// We made improvements on Brent's method, adapted from the Boost C++ library, for
// our purposes. First, we introduced a "guess" as the initial value, when previously
// the initial value was always set to the input "max" bracket value. We also introduced
// the option to take small gradient steps in the case that the initial algorithm fails
// to improve on the current best argmin value (turned on with use_gradients).
namespace Optimization {

// Options for optimization method.
enum class OptimizationMethod {
  BrentOptimization,
  BrentOptimizationWithGradients,
  GradientAscentOptimization,
  LogSpaceGradientAscentOptimization,
  NewtonOptimization
};

class OptimizationMethodEnum
    : public EnumWrapper<OptimizationMethod, size_t, 5,
                         OptimizationMethod::BrentOptimization,
                         OptimizationMethod::NewtonOptimization> {
 public:
  static inline const std::string Prefix = "OptimizationMethod";
  static inline const Array<std::string> Labels = {
      {"BrentOptimization", "BrentOptimizationWithGradients",
       "GradientAscentOptimization", "LogSpaceGradientAscentOptimization",
       "NewtonOptimization"}};

  static std::string ToString(const OptimizationMethod e) {
    std::stringstream ss;
    ss << Prefix << "::" << Labels[e];
    return ss.str();
  }
  friend std::ostream &operator<<(std::ostream &os, const OptimizationMethod e) {
    os << ToString(e);
    return os;
  }
};

// Adapted from https://www.boost.org/doc/libs/1_73_0/boost/math/tools/minima.hpp
template <bool use_gradients, class F, class T>
std::tuple<T, T> BrentMinimize(F f, T guess, T min, T max, int significant_digits,
                               size_t max_iter, T step_size) {
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

  // golden ratio, don't need too much precision here!
  static const T golden = 0.3819660f;

  w = v = x = guess;
  if constexpr (use_gradients) {
    fw = fv = fx = f(x).first;
  } else {
    fw = fv = fx = f(x);
  }
  delta2 = delta = 0;

  size_t count = max_iter;

  do {
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
          delta = (mid - x) < 0 ? static_cast<T>(-fabs(fract1))
                                : static_cast<T>(fabs(fract1));
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
    if constexpr (use_gradients) {
      fu = f(u).first;
    } else {
      fu = f(u);
    }

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
    } else if constexpr (use_gradients) {
      // Considering update using gradient descent
      T f_prime_x = f(x).second;
      T u_ = x - step_size * f_prime_x;
      T fu_ = f(u_).first;
      if (fu_ <= fx) {
        // good new point using gradient is an improvement!
        // update brackets (previous guess becomes the new outer bracket):
        if (u_ >= x)
          min = x;
        else
          max = x;
        // update control points:
        v = w;
        w = x;
        x = u_;
        fv = fw;
        fw = fx;
        fx = fu_;
      }
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
  } while (--count);  // countdown until max iterations.

  max_iter -= count;

  return std::make_tuple(x, fx);
}

double GradientAscent(std::function<DoublePair(double)> f_and_f_prime, double x,
                      const int significant_digits, const double step_size,
                      const double min_x, const size_t max_iter);

double LogSpaceGradientAscent(std::function<DoublePair(double)> f_and_f_prime, double x,
                              const int significant_digits,
                              const double log_space_step_size, const double min_x,
                              const size_t max_iter);

double NewtonRaphsonOptimization(
    std::function<std::tuple<double, double, double>(double)> f_and_derivatives,
    double x, const int significant_digits, const double epsilon, const double min_x,
    const double max_x, const size_t max_iter);
}  // namespace Optimization
