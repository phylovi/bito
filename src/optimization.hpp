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
std::ostream &myout = std::cout;
// std::ostream myout = std::ostream("/dev/null");

// #include <boost/iostreams/stream.hpp>
// #include <boost/iostreams/device/null.hpp>
using Func = std::function<double(double)>;
using FuncAndOneDerivative = std::function<DoublePair(double)>;
using FuncAndTwoDerivatives = std::function<std::tuple<double, double, double>(double)>;

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

// TODO: REMOVE THIS BEFORE MERGING WITH MAIN.
// Optimization tool global benchmarks.
int iterations = 0;
bool is_converged = false;

// Copied from https://www.boost.org/doc/libs/1_73_0/boost/math/tools/minima.hpp
template <class F, class T>
std::tuple<T, T, std::pair<std::vector<T>, std::vector<T>>> BrentMinimize(F f, T min, T max,
                                                             int significant_digits,
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

  // golden ratio, don't need too much precision here!
  static const T golden = 0.3819660f;

  x = w = v = max;
  fw = fv = fx = f(x);
  delta2 = delta = 0;

  size_t count = max_iter;

  std::vector<T> path_x;
  std::vector<T> path_fx;

  path_x.push_back(x);
  path_fx.push_back(fx);

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

    iterations++;
  } while (--count);  // countdown until max iterations.

  max_iter -= count;

  return std::make_tuple(x, fx, std::make_pair(path_x, path_fx));
}

// TODO: REMOVE THIS BEFORE MERGING WITH MAIN.
// Duplicate version with extra notes by drich.
// Copied from https://www.boost.org/doc/libs/1_73_0/boost/math/tools/minima.hpp
template <class F, class T>
std::pair<T, T> BrentMinimize_FromMe(F f, T min, T max, int significant_digits,
                                     size_t max_iter) {
  T tolerance = static_cast<T>(ldexp(1.0, 1 - significant_digits));
  T x;                 // minima so far
  T w;                 // second best point
  T v;                 // previous value of w
  T u;                 // most recent evaluation point
  T delta;             // The distance moved in the last step
  T delta2;            // The distance moved in the step before last
  T fu, fv, fw, fx;    // function evaluations at u, v, w, x
  T mid;               // midpoint of min and max
  T fract1, fract2;    // minimal relative movement in x
  T abs_min, abs_max;  // absolute min and min allowed.

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
      if ((f(max) - f(mid) >= -brackets_tol) == false) {
        max = (max + abs_max) / 2;
        max_updated = true;
      };
      // update min bracket if it does not contain optimal value
      if ((f(min) - f(mid) >= -brackets_tol) == false) {
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
      brackets_contain_optimal =
          ((f(max) - f(mid) >= -brackets_tol) && (f(min) - f(mid) >= -brackets_tol));
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

// Boost version of the TOM748 Rootfinding algorithm.
// Copied from https://www.boost.org/doc/libs/1_73_0/boost/math/tools/minima.hpp
// template <class F, class T>
namespace TOMS_748 {
// Returns -1, 1 or 0, depending on the sign.
template <class T>
inline int Sign(T &val) {
  if (val == 0) {
    return 0;
  }
  // signbit return true if val is negative.
  return (std::signbit(val) == false) ? 1 : -1;
}

// Checks if distance between brackets are less than tolerance.
template <class T>
inline bool TestTol(T &min, T &max, T &tol) {
  return (max - min) < tol;
};

// Return True if values are closer in distance than tolerance.
template <class T>
inline bool TestNotDistinct(T &a, T &b) {
  T min_diff = std::numeric_limits<T>::min() * 32;
  return (fabs(a - b) < min_diff);
};

// Return True if not all values are distinct.
template <class T>
inline bool TestNotAllDistinct(T &a, T &b, T &c, T &d, T &e) {
  T min_diff = std::numeric_limits<T>::min() * 32;
  bool not_all_distinct =
      (TestNotDistinct(a, b) || TestNotDistinct(a, d) || TestNotDistinct(a, e) ||
       TestNotDistinct(b, d) || TestNotDistinct(b, e) || TestNotDistinct(d, e));
  return not_all_distinct;
};

// Perform division if no over/underflow error occurs.
// Otherwise, return r.
template <class T>
inline T SafeDivision(T num, T denom, T r) {
  //
  // return num / denom without overflow,
  // return r if overflow would occur.
  //
  const T max_value = std::numeric_limits<T>::max();
  if (fabs(denom) < 1) {
    if (fabs(denom * max_value) <= fabs(num)) {
      return r;
    }
  }
  return num / denom;
}

// Update interval brackets [a,b] with point c, to either [a,c] or [c,b].
// Results stored in place.
// d and fd updated to the point removed from interval.
template <class F, class T>
inline void Bracket(F f, T &a, T &b, T &c, T &fa, T &fb, T &d, T &fd) {
  //
  // Given a point c inside the existing enclosing interval
  // [a, b] sets a = c if f(c) == 0, otherwise finds the new
  // enclosing interval: either [a, c] or [c, b] and sets
  // d and fd to the point that has just been removed from
  // the interval.  In other words d is the third best guess
  // to the root.
  //
  T tol = std::numeric_limits<T>::epsilon() * 2;
  //
  // If the interval [a,b] is very small, or if c is too close
  // to one end of the interval then we need to adjust the
  // location of c accordingly:
  //
  if ((b - a) < 2 * tol * a) {
    c = a + (b - a) / 2;
  } else if (c <= a + fabs(a) * tol) {
    c = a + fabs(a) * tol;
  } else if (c >= b - fabs(b) * tol) {
    c = b - fabs(a) * tol;
  }
  T fc = f(c);
  // If f(c) is zero, then we have found the exact root.
  if (fc == 0) {
    a = c;
    fa = 0;
    d = 0;
    fd = 0;
    return;
  }
  // Otherwise, update the interval to which side of c contains the root.
  if (Sign(fa) * Sign(fc) < 0) {
    d = b;
    fd = fb;
    b = c;
    fb = fc;
  } else {
    d = a;
    fd = fa;
    a = c;
    fa = fc;
  }
};

// Performs step of the secant method.
template <class T>
inline T SecantInterpolate(const T &a, const T &b, const T &fa, const T &fb) {
  //
  // Performs standard secant interpolation of [a,b] given
  // function evaluations f(a) and f(b).  Performs a bisection
  // if secant interpolation would leave us very close to either
  // a or b.  Rationale: we only call this function when at least
  // one other form of interpolation has already failed, so we know
  // that the function is unlikely to be smooth with a root very
  // close to a or b.
  //
  T tol = std::numeric_limits<T>::epsilon() * 5;
  T c = a - (fa / (fb - fa)) * (b - a);
  if ((c <= a + fabs(a) * tol) || (c >= b - fabs(b) * tol)) {
    return (a + b) / 2;
  }
  return c;
};

// Performs step of the quadratic interpolation method.
template <class T>
inline T QuadtraticInterpolate(const T &a, const T &b, const T &d, const T &fa,
                               const T &fb, const T &fd, size_t count) {
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
  T max_value = std::numeric_limits<T>::max();
  T zero = T(0);

  T qA = fd - fb;
  T dA = d - b;
  T A = SafeDivision(qA, dA, max_value);
  T qB = fb - fa;
  T dB = b - a;
  T B = SafeDivision(qB, dB, max_value);
  A = SafeDivision(T(A - B), T(d - a), zero);

  if (a == 0) {
    // failure to determine coefficients, try a secant step:
    return SecantInterpolate(a, b, fa, fb);
  }
  //
  // Determine the starting point of the Newton steps:
  //
  T c;
  if (Sign(A) * Sign(fa) > 0) {
    c = a;
  } else {
    c = b;
  }
  //
  // Take the Newton steps:
  //
  for (unsigned i = 1; i <= count; ++i) {
    // c -= SafeDivision(B * c, (B + A * (2 * c - a - b)), 1 + c - a);
    c -= SafeDivision(T(fa + (B + A * (c - b)) * (c - a)), T(B + A * (2 * c - a - b)),
                      T(1 + c - a));
  }
  if ((c <= a) || (c >= b)) {
    // Oops, failure, try a secant step:
    c = SecantInterpolate(a, b, fa, fb);
  }
  return c;
};

// Performs step of the cubic interpolation method.
template <class T>
inline T CubicInterpolate(const T &a, const T &b, const T &d, const T &e, const T &fa,
                          const T &fb, const T &fd, const T &fe) {
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
  T q11 = (d - e) * fd / (fe - fd);
  T q21 = (b - d) * fb / (fd - fb);
  T q31 = (a - b) * fa / (fb - fa);
  T d21 = (b - d) * fd / (fd - fb);
  T d31 = (a - b) * fb / (fb - fa);

  T q22 = (d21 - q11) * fb / (fe - fb);
  T q32 = (d31 - q21) * fa / (fd - fa);
  T d32 = (d31 - q21) * fd / (fd - fa);
  T q33 = (d32 - q22) * fa / (fe - fa);
  T c = q31 + q32 + q33 + a;

  // If step goes out of bounds of brackets, fall back to quadratic interpolation:
  if ((c <= a) || (c >= b)) {
    c = QuadtraticInterpolate(a, b, d, fa, fb, fd, 3);
  }

  return c;
};

// Performs step of the bijection method.
template <class T>
inline void Bijection(const T &a, const T &b, const T &d, const T &e, const T &fa,
                      const T &fb, const T &fd, const T &fe) {
  e = d;
  fe = fd;
  T mid = a + ((b - a) / 2);
  return mid;
};

// Main method.
// Return is a bracketed range containing the root.
template <class F, class T>
std::pair<T, T> RootFinder(F f, const T &min, const T &max, int &significant_digits,
                           size_t &max_iter) {
  //
  // Main entry point and logic for Toms Algorithm 748 root finder.
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
  const T mu = 0.5f;
  // Iteration countdown counter
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
  if (TestTol(a, b, tol) || (fa == 0) || (fb == 0)) {
    max_iter = 0;
    if (fa == 0) {
      b = a;
    } else if (fb == 0) {
      a = b;
    }
    // Meta data
    iterations = 0;
    is_converged = true;
    // Return ordered (min, max) span that contain the root.
    return std::make_pair(a, b);
  }

  // If f(a) and f(b) properly bracket root, then they should have opposite signs.
  std::cout << "TOMS_748:: f[" << a << "," << b << "] = (" << fa << "," << fb << ")"
            << std::endl;
  Assert(Sign(fa) * Sign(fb) < 0, "TOMS_748: <min> and <max> don't bracket the root.");

  // Initialize dummy vals for fd, e, fe.
  fe = e = fd = 1e5f;

  // On first step, use secant method.
  if (fa != 0) {
    SecantInterpolate(a, b, fa, fb);
    Bracket(f, a, b, c, fa, fb, d, fd);
    count--;
  }
  // On second step, use quadratic interpolation method.
  if ((fa != 0) && (count != 0) && !TestTol(a, b, tol)) {
    c = QuadtraticInterpolate(a, b, d, fa, fb, fd, 2);
    e = d;
    fe = fd;
    Bracket(f, a, b, c, fa, fb, d, fd);
    count--;
  }

  // Loop until convergence.
  while (((count == 0) || (fa == 0) || TestTol(a, b, tol)) == false) {
    // Save brackets.
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
    // (1) Take an iterpolation step.
    // Check if fa,fb,fc,fd are all distinct (greater than min_diff apart).
    bool prof = TestNotAllDistinct(fa, fb, fc, fd, fe);
    // If not all values are distinct, then fall back to quadratic interpolation.
    if (prof) {
      c = QuadtraticInterpolate(a, b, d, fa, fb, fd, 2);
    }
    // Otherwise, use cubic interpolation.
    else {
      c = CubicInterpolate(a, b, d, e, fa, fb, fd, fe);
    }
    // Re-bracket and check for convergence.
    Bracket(f, a, b, c, fa, fb, d, fd);
    count--;
    if ((count == 0) || (fa == 0) || TestTol(a, b, tol)) {
      break;
    }
    // (2) Take another interpolation step.
    // Check if fa,fb,fc,fd are all distinct (greater than min_diff apart).
    prof = TestNotAllDistinct(fa, fb, fc, fd, fe);
    // If not all values are distinct, then fall back to quadratic interpolation.
    if (prof) {
      c = QuadtraticInterpolate(a, b, d, fa, fb, fd, 3);
    }
    // Otherwise, use cubic interpolation.
    else {
      c = CubicInterpolate(a, b, d, e, fa, fb, fd, fe);
    }
    // Re-bracket and check for convergence.
    Bracket(f, a, b, c, fa, fb, d, fd);
    count--;
    if ((count == 0) || (fa == 0) || TestTol(a, b, tol)) {
      break;
    }
    // (3) Take a double-length secant step.
    if (fabs(fa) < fabs(fb)) {
      u = a;
      fu = fa;
    } else {
      u = b;
      fu = fb;
    }
    c = u - 2 * (fu / (fb - fa)) * (b - a);
    if (fabs(c - u) > (b - a) / 2) {
      c = a + (b - a) / 2;
    }
    // Re-bracket and check for convergence.
    Bracket(f, a, b, c, fa, fb, d, fd);
    count--;
    if ((count == 0) || (fa == 0) || TestTol(a, b, tol)) {
      break;
    }
    // (4) If not converging quickly enough, then take a bijection step.
    // If search space hasn't shrunk by at least a factor of mu, then use bijection.
    if ((b - a) >= mu * (b0 - a0)) {
      continue;
    }
    // Re-bracket again on a bisection.
    e = d;
    fe = fd;
    T mid = a + ((b - a) / 2);
    Bracket(f, a, b, mid, fa, fb, d, fd);
    count--;
  };

  // If exact root was found, collapse to single point.
  if (fa == 0) {
    b = a;
  } else if (fb == 0) {
    a = b;
  }

  // Meta results
  iterations = max_iter - count;
  is_converged = TestTol(a, b, tol);
  // Total iterations taken.
  max_iter -= count;
  // Return the span containing the root.
  return std::make_pair(a, b);
};

// Helper function to convert minimization function to root-finding function.
template <class F, class T>
std::pair<T, T> Minimize(F f, std::function<std::pair<T, T>(T)> f_and_df, const T &min,
                         const T &max, int &significant_digits, size_t &max_iter) {
  // Treat derivative function as primary function.
  auto df = [&f_and_df](T x) {
    auto [fx, dfx] = f_and_df(x);
    return dfx;
  };

  // Returns a span containing optimal x value.
  auto [results_min, results_max] =
      RootFinder<F, T>(df, min, max, significant_digits, max_iter);
  auto result = (results_min + ((results_max - results_min) / 2));

  // return average and function evaluated at average.
  return std::make_pair(result, f(result));
};

};  // namespace TOMS_748

//
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
std::tuple<double, double, std::pair<std::vector<double>, std::vector<double>>>
NewtonRaphsonOptimization(
    std::function<std::tuple<double, double, double>(double)> f_and_derivatives,
    double x, const double tolerance, const double epsilon, const double min_x,
    const double max_x, const size_t max_iter) {
  size_t iter_idx = 0;
  double new_x, delta;
  double damp_const = 0.005;

  std::vector<double> path_x;
  std::vector<double> path_fx;

  while (true) {
    auto [f_x, f_prime_x, f_double_prime_x] = f_and_derivatives(x);
    path_x.push_back(x);
    path_fx.push_back(f_x);

    if (fabs(f_double_prime_x) < epsilon) {
      return std::make_tuple(x, f_x, std::make_pair(path_x, path_fx));
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
      return std::make_tuple(x, f_x, std::make_pair(path_x, path_fx));
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

// All below adapted from boost
//
template <class T>
inline int sign(const T z) {
  return (z == 0) ? 0 : std::signbit(z) ? -1 : 1;
}

template <class F, class T>
void handle_zero_derivative(F f, T last_f1, const T f1, T delta, T result, T guess,
                            const T min, const T max) {
  if (last_f1 == 0) {
    // this must be the first iteration, pretend that we had a
    // previous one at either min or max:
    if (result == min) {
      guess = max;
    } else {
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

// TODO: REMOVE THIS BEFORE MERGING WITH MAIN.
// This just selects a various minimization functions for comparison.
template <class F, class T>
std::pair<T, T> Optimize_RunAll(F f, F neg_f, FuncAndOneDerivative f_and_df,
                                FuncAndTwoDerivatives f_and_df_and_ddf, T min, T max,
                                T log_min, T log_max, int significant_digits,
                                size_t max_iter) {
  int iters;
  bool converged;
  std::pair<T, T> results;
  std::pair<T, T> results_final;
  T abs_min = std::numeric_limits<T>::min();

  auto report_results = [&](std::string method_name, bool use_log_lengths = false,
                            bool use_results = false) {
    iters = iterations;
    converged = is_converged;

    auto x = (use_log_lengths ? exp(results.first) : results.first);
    auto fx = (use_log_lengths ? results.second : results.second);

    myout << method_name << ": (x,fx) = [" << x << ", " << fx << "]" << std::endl;
    myout << method_name << ": (x-,x+) = [" << f(abs_min) << ", " << f(2 * x) << "]"
          << std::endl;
    myout << method_name << ": "
          << "iters: " << iterations << ", converged?: " << is_converged << std::endl;
  };

  auto df = [&f_and_df](T x) {
    auto [fx, dfx] = f_and_df(x);
    return dfx;
  };

  // Print out ordered pair to get shape of function.
  auto print_func_pairs = [&](F func, std::string function_name,
                              bool use_log_lengths = false, int num_steps = 100) {
    T my_min = (use_log_lengths ? log_min : min);
    T my_abs_min = (use_log_lengths ? log(abs_min) : abs_min);
    my_min = my_abs_min;
    T my_max = (use_log_lengths ? log_max : max);
    T step_size = (max - min) / T(num_steps);
    std::cout << function_name << ":: (Min,Max): " << my_min << ", " << my_max
              << ", Dist: " << my_max - my_min << ", Stepsize: " << step_size
              << std::endl;
    int per_line = 5;
    int step_cnt = 1;
    for (T i = my_min; i < my_max; i += step_size) {
      std::cout << "[" << step_cnt << "](" << i << ", " << func(i) << ") ";

      if (step_cnt % per_line == 0) {
        std::cout << std::endl;
      }
      step_cnt++;
      if (step_cnt > 100) {
        break;
      }
    }
    std::cout << std::endl;
  };

  std::cout << std::endl << std::endl;
  std::cout << "OPTIMIZE_RUNALL:: f[" << min << ", " << max << "] = (" << f(min) << ", "
            << f(max) << ")" << std::endl;
  std::cout << "OPTIMIZE_RUNALL:: df[" << min << ", " << max << "] = (" << df(min)
            << ", " << df(max) << ")" << std::endl;
  std::cout << "OPTIMIZE_RUNALL:: df[" << abs_min << "] = (" << df(abs_min) << ")"
            << std::endl;

  print_func_pairs(neg_f, "neg. log_likelihood w/ log_lengths", true);
  print_func_pairs(f, "log_likelihood");
  print_func_pairs(df, "deriv. log_likelihood");

  // Run Brent Minimization.
  // results = BrentMinimize<F, T>(neg_f, log_min, log_max, significant_digits,
  // max_iter); report_results("Brent_Minimize", true, true);

  // Run TOMS748 Minimization.
  results =
      TOMS_748::Minimize<F, T>(f, f_and_df, min, max, significant_digits, max_iter);
  report_results("TOMS_748", false);

  // Run Newton-Raphson Optimization.
  // results = NewtonRaphsonOptimization();
  // report_results("NewtonRaphson", false);

  return results;
}

}  // namespace Optimization
