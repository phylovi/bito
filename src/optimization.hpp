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
          delta = (mid - x) < 0 ? (T)-fabs(fract1) : (T)fabs(fract1);
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

DoublePair NewtonRaphsonOptimization(
    std::function<std::tuple<double, double, double>(double)> f_and_derivatives,
    double x, const double tolerance, const double epsilon, const double min_x,
    const double max_x, const size_t max_iter) {
  size_t iter_idx = 0;
  double new_x, delta;

  while (true) {
    auto [f_x, f_prime_x, f_double_prime_x] = f_and_derivatives(x);
    new_x = x - f_prime_x / f_double_prime_x;
    new_x = fmax(new_x, min_x);
    if (new_x >= max_x) {
      new_x = x - 0.5 * (x - max_x);
    }

    delta = fabs(x - new_x);

    if (delta < tolerance || fabs(f_double_prime_x) < epsilon || iter_idx == max_iter) {
      return {x, f_x};
    } else {
      x = new_x;
    }
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

// All below adapted from boost
//
template <class T>
inline int sign(const T z) {
  return (z == 0) ? 0 : std::signbit(z) ? -1 : 1;
}

template <class F, class T>
void handle_zero_derivative(F f, T last_f0, const T f0, T delta, T result, T guess,
                            const T min, const T max) {
  // if(last_f0 == 0)
  // {
  //    // this must be the first iteration, pretend that we had a
  //    // previous one at either min or max:
  //    if(result == min)
  //    {
  //       guess = max;
  //    }
  //    else
  //    {
  //       guess = min;
  //    }
  //    unpack_0(f(guess), last_f0);
  //    delta = guess - result;
  // }
  if (sign(last_f0) * sign(f0) < 0) {
    // we've crossed over so move in opposite direction to last step:
    if (delta < 0) {
      delta = (result - min) / 2;
    } else {
      delta = (result - max) / 2;
    }
  } else {
    // move in same direction as last step:
    if (delta < 0) {
      delta = (result - max) / 2;
    } else {
      delta = (result - min) / 2;
    }
  }
}

template <class F, class T>
std::pair<T, T> NewtonRaphsonIterate(F f, T guess, T min, T max, int significant_digits,
                                     size_t max_iter) {
  T f1(0), f2, last_f1(0);
  T result = guess;

  T factor = static_cast<T>(ldexp(1.0, 1 - significant_digits));
  T delta = 1;
  T delta1 = std::numeric_limits<double>::max();
  T delta2 = std::numeric_limits<double>::max();

  size_t count = max_iter;

  do {
    last_f1 = f1;
    delta2 = delta1;
    delta1 = delta;
    auto [f0, f1, f2] = f(result);
    if (0 == f1) break;
    if (f2 == 0) {
      // Oops zero derivative!!!
      handle_zero_derivative(f, last_f1, f1, delta, result, guess, min, max);
    } else {
      delta = f1 / f2;
    }
    if (fabs(delta * 2) > fabs(delta2)) {
      // last two steps haven't converged, try bisection:
      delta = (delta > 0) ? (result - min) / 2 : (result - max) / 2;
    }
    guess = result;
    result -= delta;
    if (result <= min) {
      delta = 0.5F * (guess - min);
      result = guess - delta;
      if ((result == min) || (result == max)) break;
    } else if (result >= max) {
      delta = 0.5F * (guess - max);
      result = guess - delta;
      if ((result == min) || (result == max)) break;
    }
    // update brackets:
    if (delta > 0)
      max = guess;
    else
      min = guess;
  } while (--count && (fabs(result * factor) < fabs(delta)));

  max_iter -= count;

  return std::make_pair(result, std::get<0>(f(result)));
}
}  // namespace Optimization
