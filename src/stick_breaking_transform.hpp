// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_STICK_BREAKING_TRANSFORM_HPP_
#define SRC_STICK_BREAKING_TRANSFORM_HPP_

#include "eigen_sugar.hpp"

class Transform {
 public:
  virtual ~Transform() = default;

  virtual EigenVectorXd operator()(EigenVectorXd const& x) const = 0;

  virtual EigenVectorXd inverse(EigenVectorXd const& y) const = 0;

  virtual double log_abs_det_jacobian(const EigenVectorXd& x,
                                      const EigenVectorXd& y) const = 0;
};

class StickBreakingTransform : public Transform {
  // The stick breaking procedure as defined in Stan
  // https://mc-stan.org/docs/2_26/reference-manual/simplex-transform-section.html
 public:
  virtual ~StickBreakingTransform() = default;

  EigenVectorXd operator()(EigenVectorXd const& x) const;

  EigenVectorXd inverse(EigenVectorXd const& y) const;

  double log_abs_det_jacobian(const EigenVectorXd& x, const EigenVectorXd& y) const;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("BreakingStickTransform") {
  StickBreakingTransform a;
  EigenVectorXd y(3);
  y << 1., 2., 3.;
  EigenVectorXd x_expected(3);
  // x_expected =
  // torch.distributions.StickBreakingTransform()(torch.tensor([1., 2., 3.]))
  x_expected << 0.475367, 0.412879, 0.106454, 0.00530004;
  EigenVectorXd x = a(y);
  CheckVectorXdEquality(x, x_expected, 1.e-5);
  EigenVectorXd yy = a.inverse(x);
  CheckVectorXdEquality(y, yy, 1e-5);
  // log_abs_det_jacobian_expected =
  // torch.distributions.StickBreakingTransform().log_abs_det_jacobian(y,x)
  double log_abs_det_jacobian_expected = -9.108352;
  CHECK_EQ(a.log_abs_det_jacobian(x, y), log_abs_det_jacobian_expected, 1.e-5);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_STICK_BREAKING_TRANSFORM_HPP_
