// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Put Eigen "common" code here.

#ifndef SRC_EIGEN_SUGAR_HPP_
#define SRC_EIGEN_SUGAR_HPP_

#include <Eigen/Dense>
#include <fstream>

#include "sugar.hpp"

using EigenVectorXd = Eigen::VectorXd;
using EigenVectorXi = Eigen::VectorXi;
using EigenMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using EigenMatrixXb = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
using EigenVectorXdRef = Eigen::Ref<EigenVectorXd>;
using EigenMatrixXdRef = Eigen::Ref<EigenMatrixXd>;
using EigenConstVectorXdRef = Eigen::Ref<const EigenVectorXd>;
using EigenConstMatrixXdRef = Eigen::Ref<const EigenMatrixXd>;
using EigenArrayXb = Eigen::Array<bool, Eigen::Dynamic, 1>;
using EigenArrayXbRef = Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, 1>>;

const static Eigen::IOFormat EigenCSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols,
                                            ", ", "\n");
// Write an Eigen object to a CSV file.
template <class EigenType>
void EigenToCSV(const std::string &file_path, EigenType eigen_object) {
  std::ofstream file(file_path.c_str());
  file << eigen_object.format(EigenCSVFormat) << std::endl;
  if (file.bad()) {
    Failwith("Failure writing to " + file_path);
  }
}

// Convert each entry of a std::vector<T> to double using a function f and store in an
// EigenVectorXd.
template <typename T>
EigenVectorXd EigenVectorXdOfStdVectorT(const std::vector<T> &v,
                                        const std::function<double(const T &)> &f) {
  EigenVectorXd results(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    results[i] = f(v[i]);
  }
  return results;
}

#ifdef DOCTEST_LIBRARY_INCLUDED
void CheckVectorXdEquality(double value, const EigenVectorXd v, double tolerance) {
  for (size_t i = 0; i < v.size(); i++) {
    CHECK_LT(fabs(value - v[i]), tolerance);
  }
};

void CheckVectorXdEquality(const EigenVectorXd v1, const EigenVectorXd v2,
                           double tolerance) {
  CHECK_EQ(v1.size(), v2.size());
  for (size_t i = 0; i < v1.size(); i++) {
    CHECK_LT(fabs(v1[i] - v2[i]), tolerance);
  }
};

void CheckVectorXdEqualityAfterSorting(const EigenVectorXdRef v1,
                                       const EigenVectorXdRef v2, double tolerance) {
  EigenVectorXd v1_sorted = v1;
  EigenVectorXd v2_sorted = v2;
  std::sort(v1_sorted.begin(), v1_sorted.end());
  std::sort(v2_sorted.begin(), v2_sorted.end());
  CheckVectorXdEquality(v1_sorted, v2_sorted, tolerance);
};

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_EIGEN_SUGAR_HPP_
