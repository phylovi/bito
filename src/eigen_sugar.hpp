// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Put Eigen "common" code here.

#ifndef SRC_EIGEN_SUGAR_HPP_
#define SRC_EIGEN_SUGAR_HPP_

#include <Eigen/Dense>

using EigenMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using EigenVectorXdRef = Eigen::Ref<Eigen::VectorXd>;
using EigenMatrixXdRef = Eigen::Ref<EigenMatrixXd>;

#ifdef DOCTEST_LIBRARY_INCLUDED
void CheckVectorXdEquality(double value, const Eigen::VectorXd v,
                           double tolerance) {
  for (size_t i = 0; i < v.size(); i++) {
    CHECK_LT(fabs(value - v[i]), 0.0001);
  }
};

void CheckVectorXdEquality(const Eigen::VectorXd v1, const Eigen::VectorXd v2,
                           double tolerance) {
  for (size_t i = 0; i < v1.size(); i++) {
    CHECK_LT(fabs(v1[i] - v2[i]), 0.0001);
  }
};
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_EIGEN_SUGAR_HPP_
