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

#endif  // SRC_EIGEN_SUGAR_HPP_
