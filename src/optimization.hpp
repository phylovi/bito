// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <cmath>
#include <functional>
#include <numeric>
#include <utility>

#include "sugar.hpp"

namespace Optimization {
  // test stream
  std::ofstream myout = std::cout;
  // std::ofstream myout = std::ofstream("/dev/null");

// #include <boost/iostreams/stream.hpp>
// #include <boost/iostreams/device/null.hpp>

// TODO: template using?
using FuncAndOneDerivative = std::function<DoublePair(double)>;
using FuncAndTwoDerivatives = std::function<std::tuple<double, double, double>(double)>;

// Various Minimizing Functions.
template <class F, class T>
std::pair<T, T> BrentMinimize_FromBoost(F f, T min, T max, int significant_digits,
                                      size_t max_iter);
template <class F, class T>
std::pair<T, T> BrentMinimize_FromMe(F f, T min, T max, int significant_digits,
                              size_t max_iter);   

// Write ordered pairs from function at intervals.
template <class F, class T>
void FunctionToCSV(F f, T min, T max, T step_size, std::string filepath) {
  std::ofstream os = std::ofstream(filepath);
  os << "# x,f(x)" << std::endl;
  for (int x = min; x < max; x += step_size) {
    os << x << "," << f(x) << std::endl;
  }
  os.close();
};

// Optimization tool benchmarks.
int iterations = 0;
bool is_converged = false;

// TODO: temp selector / comparator
// This just selects a various minimization functions for comparison.
template <class F, class T>
std::pair<T, T> Minimize(F f, T min, T max, int significant_digits,
                         size_t max_iter) {

  // Brent from Boost.
  int boost_iter;
  bool boost_converged;
  auto results_boost 
      = BrentMinimize_FromBoost<F,T>(f, min, max, significant_digits, max_iter);
  boost_iter = iterations;
  boost_converged = is_converged;

  // Brent from me.
  int my_iter;
  bool my_converged;
  auto results_me 
      = BrentMinimize_FromMe<F,T>(f, min, max, significant_digits, max_iter);
  my_iter = iterations;
  my_converged = is_converged;

  // Output comparison.
  myout << "BRENT_BOOST: " << "(f, fx) = [" << results_boost.first << "," << results_boost.second << "]" << std::endl;
  myout << "BRENT_BOOST: " << "iters: " << boost_iter << ", converged?: " << is_converged << std::endl;

  myout << "BRENT_ME: " << "(f, fx) = [" << results_me.first << "," << results_me.second << "]" << std::endl;
  myout << "BRENT_ME: " << "iters: " << my_iter << ", converged?: " << is_converged << std::endl;

  FunctionToCSV<F,T>(f, min, max, (max - min)/10000, "_ignore/function_shape_from_optimization.cs");

  return results_boost;
}

// Copied from https://www.boost.org/doc/libs/1_73_0/boost/math/tools/minima.hpp
template <class F, class T>
std::pair<T, T> BrentMinimize_FromBoost(F f, T min, T max, int significant_digits,
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

  iterations = 0;
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
    bool use_bisection = true; // only triggers if IQI fails.
    if (fabs(delta2) > fract1) {
      // try and construct a parabolic fit:
      // TODO: where is secant method? what happens when x=v=w and fx=fv=fw?
      T r = (x - w) * (fx - fv);
      T q = (x - v) * (fx - fw);
      T p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) p = -p;
      q = fabs(q);
      T td = delta2;
      delta2 = delta;
      // determine whether a parabolic step is acceptable or not:
      // must fail all three threshold tests to be accepted.
      if ( ((fabs(p) >= fabs(q * td / 2)) == false) 
           && ((p <= q * (min - x)) == false)
           && ((p >= q * (max - x)) == false) ) {
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

    // ** Golden Bisection Method (this is an optimization of traditional Bisection Method)
    if (use_bisection) {
      // golden section:
      delta2 = (x >= mid) ? min - x : max - x;
      delta = golden * delta2;
    }

    // ** Update current position:
    u = (fabs(delta) >= fract1) ? 
        T(x + delta) : (delta > 0 ? T(x + fabs(fract1)) : T(x - fabs(fract1)));
    fu = f(u);
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
    } 
    else {
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
      } 
      else if ((fu <= fv) || (v == x) || (v == w)) {
        // third best:
        v = u;
        fv = fu;
      }
    }

  iterations++;
  } while (--count); // countdown until max iterations.

  // What does this for?
  max_iter -= count;

  return std::make_pair(x, fx);
}

// Copied from https://www.boost.org/doc/libs/1_73_0/boost/math/tools/minima.hpp
template <class F, class T>
std::pair<T, T> BrentMinimize_FromMe(F f, T min, T max, int significant_digits,
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
  T abs_min, abs_max; // absolute min and min allowed.

  static const T golden =
      0.3819660f;  // golden ratio, don't need too much precision here!

  x = w = v = max;
  fw = fv = fx = f(x);
  delta2 = delta = 0;

  size_t count = max_iter;

  // initialize midpoint.
  abs_min = min;
  abs_max = max;
  mid = (min + max) / 2;
  myout << "MID_INIT(" << min << "," << max << "): " << mid << std::endl;
  bool first_pass = true;
  
  do {
    // get midpoint
    mid = (min + max) / 2;

    // update brackets: expansion search via binary search.
    size_t brackets_iter = 0;
    size_t brackets_iter_max = 100;
    T brackets_tol = tolerance;
    bool brackets_contain_optimal;
    bool min_updated, max_updated;
    do {
      // update max bracket if it does not contain optimal value
      if ( (f(max) - f(mid) >= -brackets_tol) == false ) { 
        max = (max + abs_max) / 2; 
        max_updated = true;
      };
      // update min bracket if it does not contain optimal value
      if ( (f(min) - f(mid) >= -brackets_tol) == false ) { 
        min = (min + abs_min) / 2; 
        min_updated = true;
      };
      brackets_iter++;
      myout << "ITER: " << brackets_iter << " | ";
      myout << "MID: " << mid << " " << f(mid) << " | "; 
      myout << (min_updated ? "*" : "") << "MIN: " << min << " " << f(min) << " | "; 
      myout << (max_updated ? "*" : "") << "MAX: " << max << " " << f(max) << std::endl; 
      if (brackets_iter >= brackets_iter_max) { 
        myout << "SEARCH TIMED OUT!" << std::endl;
        break;
      }
      min_updated = false;
      max_updated = false;
      // continue until brackets contain optimal value (within tolerance)
      brackets_contain_optimal = ((f(max) - f(mid) >= -brackets_tol) && (f(min) - f(mid) >= -brackets_tol));
    } while (brackets_contain_optimal == false);

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

  myout << "MINIMUM FOUND: " << x << " " << fx << std::endl;
  return std::make_pair(x, fx);
}


template <class F, class T>
class TOMS_748 {
 public:

  // Checks if distance between brackets are less than tolerance.
  inline bool TestTol(T &min, T &max, T &tol) {
    return (max - min) < tol;
  };

  // Return True if values are closer in distance than tolerance.
  inline bool TestNotDistinct(T &a, T &b, T &min_diff) {
    return fabs(fa - fb) < min_diff;
  }

  // TODO: Implement
  inline T SafeDivision(T num, T denom, T r) {
    //
    // return num / denom without overflow,
    // return r if overflow would occur.
    //
  }

  // TODO: Implement
  inline void Bracket(F f, T &a, T &b, T &d,
                      T &fa, T &fb, T &fd) {
    
  }

  // TODO: Implement
  // Performs the secant method.
  inline T SecantInterpolate(const T &a, const T &b, 
                             const T &fa, const T &fb) {
    //
    // Performs standard secant interpolation of [a,b] given
    // function evaluations f(a) and f(b).  Performs a bisection
    // if secant interpolation would leave us very close to either
    // a or b.  Rationale: we only call this function when at least
    // one other form of interpolation has already failed, so we know
    // that the function is unlikely to be smooth with a root very
    // close to a or b.
    //


  };

  // TODO: Implement
  // Performs the quadratic interpolation method.
  inline T QuadtraticInterpolate(const T &a, const T &b, const T &d, 
                                 const T &fa, const T &fb, const T &fd, 
                                 size_t count) {
    //
    // Performs quadratic interpolation to determine the next point,
    // takes count Newton steps to find the location of the
    // quadratic polynomial.
    //
    // Point d must lie outside of the interval [a,b], it is the third
    // best approximation to the root, after a and b.
    //
    // Note: this does not guarentee to find a root
    // inside [a, b], so we fall back to a secant step should
    // the result be out of range.
    //
    // Start by obtaining the coefficients of the quadratic polynomial:
    //


  };

  // TODO: Implement
  // Performs the cubic interpolation method.
  inline T CubicInterpolate(const T &a, const T &b, const T &d, const T &e,
                            const T &fa, const T &fb, const T &fd, const T &fe) {
    //
    // Uses inverse cubic interpolation of f(x) at points 
    // [a,b,d,e] to obtain an approximate root of f(x).
    // Points d and e lie outside the interval [a,b]
    // and are the third and forth best approximations
    // to the root that we have found so far.
    //
    // Note: this does not guarentee to find a root
    // inside [a, b], so we fall back to quadratic
    // interpolation in case of an erroneous result.
    //


  };
    
  // TODO: Implement
  static std::pair<T, T> RootFinder(
      std::function<T(T)> f,
      std::function<std::pair<T,T>(T)> f_and_one_derivative,
      const T &min, const T &max, 
      int significant_digits, size_t max_iter
  ) {  
    //
    // Main entry point and logic for Toms Algorithm 748
    // root finder.
    //

    // Vars
    T a, fa, a0;
    T b, fb, b0;
    T c, fc;
    T u, fu;
    T d, fd;
    T e, fe;
    // Constants
    T tol = static_cast<T>(ldexp(1.0, 1 - significant_digits));
    T min_diff = tol * 32; // TODO: change expression
    const T mu = 0.5f;
    // Iteration counter
    size_t count = max_iter;
    size_t iter = 0;

    // Initialize a, b, fa, fb.
    // (guessing these are the min/max vals)
    a = min;
    b = max;
    Assert(a < b, "TOMS_748: Args out of order. <min> cannot be greater than <max>.");
    fa = f(a);
    fb = f(b);

    // Check if we have converged on the root.
    if (TestTol(a,b,tol) || (fa == 0) || (fb == 0) ) {
      max_iter = 0;
      if (fa == 0) {
        b = a;
      }
      else if (fb == 0) {
        a = b;
      }
      // Return ordered (x, y=f(x)) pair at point of convergence.
      T x = (a + b)/2;
      return std::make_pair(x, f(x));
    }

    // If f(a) and f(b) properly bracket root, then they should have opposite signs.
    Assert(fa * fb < 0, "TOMS_748: <min> and <max> don't bracket the root.");

    // Initialize dummy vals for fd, e, fe.
    fe = e = fd = 1e5f;

    // On first step, use secant method.
    if (fa != 0) {
      SecantInterpolate(a, b, fa, fb);
      Bracket(f, a, b, c, d, fa, fb, fd);
    }
    // On second step, use quadratic interpolation method.
    if ((fa != 0) && (count != 0) && !TestTol(a,b,tol) {
      c = QuadtraticInterpolate(&a, &b, &d, &fa, &fb, &fd, 2);
    }

    while (count && (fa != 0) && !TolTest(a,b,tol)) {
      // save brackets
      a0 = a;
      b0 = b;
      //
      // Starting with the third step taken
      // we can use either quadratic or cubic interpolation.
      // Cubic interpolation requires that all four function values
      // fa, fb, fd, and fe are distinct, should that not be the case
      // then variable prof will get set to true, and we'll end up
      // taking a quadratic step instead.
      //
      T min_diff = 
      // Check if fa,fb,fc,fd all distinct.
      bool prof = (TestNotDistinct(fa, fb, min_diff) || (fabs(fa - fd) < min_diff) ||  
    }

    T dummy1, dummy2;
    return std::make_pair(dummy1, dummy2);
  };

  // Helper function to convert minimization function to root-finding function.
  template<class F, class T>
  std::pair<T, T> Minimize(
      std::function<T(T)> f,
      std::function<std::pair<T,T>(T)> f_and_one_derivative,
      std::function<std::tuple<T,T,T>(T)> f_and_two_derivatives,
      const T min, const T max,
      int significant_digits, const size_t max_iter
  ) {
    // Treat derivate function as primary function.
    auto df = [&f_and_one_derivative](T x) {
      auto [fx, dfx] = f_and_one_derivative(x);
      return dfx;
    };
    auto df_and_one_derivative = [&f_and_two_derivatives](T x) {
      auto [fx, dfx, ddfx] = f_and_two_derivatives(x);
      return std::make_pair(dfx, ddfx);
    };

    return TOMS748_RootFinder(df, df_and_one_derivative, &min, &max, significant_digits, max_iter);
  };

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
    double x, const double tolerance, const    // Main entry point and logic for Toms Algorithm 748
   // root finder.
   //double epsilon, const double min_x,
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
void handle_zero_derivative(F f, T last_f1, const T f1, T delta, T result, T guess,
                            const T min, const T max) {
  if(last_f1 == 0)
  {
      // this must be the first iteration, pretend that we had a
      // previous one at either min or max:
      if(result == min)
      {
         guess = max;
      }
      else
      {
         guess = min;
      }
      // unpack_0(f(guess), last_f0);
      last_f1 = std::get<1>(f(guess));
      delta = guess - result;
   }
  if (sign(last_f1) * sign(f1) < 0) {
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


