// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "eigen_sugar.hpp"

class Transform {
 public:
  virtual ~Transform() = default;

  virtual EigenVectorXd operator()(EigenVectorXd const& x) const = 0;

  virtual EigenVectorXd inverse(EigenVectorXd const& y) const = 0;

  virtual double log_abs_det_jacobian(const EigenVectorXd& x,
                                      const EigenVectorXd& y) const = 0;
};

class IdentityTransform : public Transform {
 public:
  virtual ~IdentityTransform() = default;

  EigenVectorXd operator()(EigenVectorXd const& x) const override { return x; };

  EigenVectorXd inverse(EigenVectorXd const& y) const override { return y; };

  double log_abs_det_jacobian(const EigenVectorXd& x,
                              const EigenVectorXd& y) const override {
    return 0;
  };
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
