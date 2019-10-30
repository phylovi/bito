// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SUBSTITUTION_MODEL_HPP_
#define SRC_SUBSTITUTION_MODEL_HPP_
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "sugar.hpp"

using EigenMatrix4d = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
using EigenVector4d = Eigen::Vector4d;

class SubstitutionModel {
 public:
  using Parameterization = std::unordered_map<std::string, Eigen::VectorXd>;

  static Eigen::VectorXd GetFromParameterization(
      Parameterization parameterization, std::string key,
      size_t expected_length) {
    auto search = parameterization.find(key);
    if (search == parameterization.end()) {
      Failwith("Model parameter " + key + " needed in model parameterization.");
    }  // else
    auto parameter = search->second;
    if (parameter.size() != expected_length) {
      Failwith("Model parameter " + key + " has length " +
               std::to_string(parameter.size()) + " but expected size was " +
               std::to_string(expected_length) + ".");
    }  // else
    return parameter;
  }

  size_t GetStateCount() const { return frequencies_.size(); }

  const std::vector<double>& GetFrequencies() { return frequencies_; }

  const std::vector<double>& GetEigenVectors() {
    UpdateEigenDecomposition();
    return evec_;
  }
  const std::vector<double>& GetInverseEigenVectors() {
    UpdateEigenDecomposition();
    return ivec_;
  }
  const std::vector<double>& GetEigenValues() {
    UpdateEigenDecomposition();
    return eval_;
  }

  virtual void SetParameters(Parameterization parameterization) = 0;

 protected:
  virtual void UpdateEigenDecomposition() = 0;

  std::vector<double> frequencies_;
  std::vector<double> evec_;
  std::vector<double> ivec_;
  std::vector<double> eval_;
};

class JCModel : public SubstitutionModel {
 public:
  JCModel() {
    frequencies_.assign(4, 0.25);
    evec_ = {1.0, 2.0, 0.0, 0.5,  1.0, -2.0, 0.5,  0.0,
             1.0, 2.0, 0.0, -0.5, 1.0, -2.0, -0.5, 0.0};

    ivec_ = {0.25, 0.25, 0.25, 0.25, 0.125, -0.125, 0.125, -0.125,
             0.0,  1.0,  0.0,  -1.0, 1.0,   0.0,    -1.0,  0.0};

    eval_ = {0.0, -1.3333333333333333, -1.3333333333333333,
             -1.3333333333333333};
  }

  void SetParameters(Parameterization parameterization) override {
    if (!parameterization.empty()) {
      Failwith("You tried to set parameters of a JC model, which has none.");
    }
  };

 protected:
  void UpdateEigenDecomposition() override {}
};

class GTRModel : public SubstitutionModel {
 public:
  GTRModel() {
    evec_.resize(16);
    ivec_.resize(16);
    eval_.resize(4);
    rates_.assign(6, 1.0);
    frequencies_.assign(4, 0.25);
  }

  void SetParameters(
      SubstitutionModel::Parameterization specification) override;

 protected:
  void UpdateEigenDecomposition() override;

 private:
  std::vector<double> rates_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <algorithm>
TEST_CASE("SubstitutionModel") {
  auto gtr_model = std::make_unique<GTRModel>();
  auto jc_model = std::make_unique<JCModel>();
  std::vector<double> evals_jc = jc_model->GetEigenValues();
  std::vector<double> evals_gtr = gtr_model->GetEigenValues();
  std::sort(evals_jc.begin(), evals_jc.end());
  std::sort(evals_gtr.begin(), evals_gtr.end());
  for (size_t i = 0; i < evals_jc.size(); i++) {
    CHECK_LT(fabs(evals_jc[i] - evals_gtr[i]), 0.0001);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SUBSTITUTION_MODEL_HPP_
