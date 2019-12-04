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
auto CheckVectorXdEquality = [](Eigen::VectorXd eval1, Eigen::VectorXd eval2,
                                double tolerance) {
  for (size_t i = 0; i < eval1.size(); i++) {
    CHECK_LT(fabs(eval1[i] - eval2[i]), 0.0001);
  }
};

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_EIGEN_SUGAR_HPP_
