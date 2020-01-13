// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Put Eigen "common" code here.

#ifndef SRC_EIGEN_SUGAR_HPP_
#define SRC_EIGEN_SUGAR_HPP_

#include <Eigen/Dense>

using EigenVectorXd = Eigen::VectorXd;
using EigenMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using EigenVectorXdRef = Eigen::Ref<EigenVectorXd>;
using EigenMatrixXdRef = Eigen::Ref<EigenMatrixXd>;
using EigenConstVectorXdRef = Eigen::Ref<const EigenVectorXd>;
using EigenConstMatrixXdRef = Eigen::Ref<const EigenMatrixXd>;

#ifdef DOCTEST_LIBRARY_INCLUDED
void CheckVectorXdEquality(double value, const EigenVectorXd v, double tolerance) {
  for (size_t i = 0; i < v.size(); i++) {
    CHECK_LT(fabs(value - v[i]), tolerance);
  }
};

void CheckVectorXdEquality(const EigenVectorXd v1, const EigenVectorXd v2,
                           double tolerance) {
  for (size_t i = 0; i < v1.size(); i++) {
    CHECK_LT(fabs(v1[i] - v2[i]), tolerance);
  }
};
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_EIGEN_SUGAR_HPP_
