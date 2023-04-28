#pragma once

#include "../src/stick_breaking_transform.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("BreakingStickTransform") {
  StickBreakingTransform a;
  EigenVectorXd y(3);
  y << 1., 2., 3.;
  EigenVectorXd x_expected(4);
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
  CHECK(a.log_abs_det_jacobian(x, y) ==
        doctest::Approx(log_abs_det_jacobian_expected).epsilon(1.e-5));
}
#endif  // DOCTEST_LIBRARY_INCLUDED
