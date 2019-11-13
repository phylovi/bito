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

using EigenMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

class SubstitutionModel {
 public:
  virtual ~SubstitutionModel() = default;

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

  const EigenMatrixXd& GetQMatrix() { return Q_; }

  const Eigen::VectorXd& GetFrequencies() { return frequencies_; }

  const EigenMatrixXd& GetEigenvectors() { return evec_; }
  const EigenMatrixXd& GetInverseEigenvectors() { return ivec_; }
  const Eigen::VectorXd& GetEigenvalues() { return eval_; }

  virtual void SetParameters(Parameterization parameterization) = 0;

  static std::unique_ptr<SubstitutionModel> OfSpecification(
      const std::string& specification);

 protected:
  Eigen::VectorXd frequencies_;
  // TODO rename these
  EigenMatrixXd evec_;
  EigenMatrixXd ivec_;
  Eigen::VectorXd eval_;
  EigenMatrixXd Q_;
};

class JC69Model : public SubstitutionModel {
 public:
  JC69Model() {
    frequencies_.resize(4);
    evec_.resize(4, 4);
    ivec_.resize(4, 4);
    eval_.resize(4);
    Q_.resize(4, 4);
    frequencies_ << 0.25, 0.25, 0.25, 0.25;
    evec_ << 1.0, 2.0, 0.0, 0.5, 1.0, -2.0, 0.5, 0.0, 1.0, 2.0, 0.0, -0.5, 1.0,
        -2.0, -0.5, 0.0;

    ivec_ << 0.25, 0.25, 0.25, 0.25, 0.125, -0.125, 0.125, -0.125, 0.0, 1.0,
        0.0, -1.0, 1.0, 0.0, -1.0, 0.0;

    eval_ << 0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333;

    Q_ << 1.0 / 3.0, -1.0, -1.0, -1.0, -1.0, 1.0 / 3.0, -1.0, -1.0, -1.0, -1.0,
        1.0 / 3.0, -1.0, -1.0, -1.0, -1.0, 1.0 / 3.0;
  }

  void SetParameters(Parameterization parameterization) override {
    if (!parameterization.empty()) {
      Failwith("You tried to set parameters of a JC model, which has none.");
    }
  }
};

class GTRModel : public SubstitutionModel {
 public:
  GTRModel() {
    // TODO can we make this in the initializer list?
    rates_.resize(6);
    frequencies_.resize(4);
    rates_ << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    frequencies_ << 0.25, 0.25, 0.25, 0.25;
    evec_.resize(4, 4);
    ivec_.resize(4, 4);
    eval_.resize(4);
    Q_.resize(4, 4);
    UpdateEigenDecomposition();
  }

  void SetParameters(SubstitutionModel::Parameterization parameters) override;

 protected:
  void UpdateEigenDecomposition();
  void UpdateQMatrix();

 private:
  Eigen::VectorXd rates_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <algorithm>
TEST_CASE("SubstitutionModel") {
  auto gtr_model = std::make_unique<GTRModel>();
  auto jc_model = std::make_unique<JC69Model>();
  Eigen::VectorXd evals_jc = jc_model->GetEigenvalues();
  Eigen::VectorXd evals_gtr = gtr_model->GetEigenvalues();
  std::sort(evals_jc.begin(), evals_jc.end());
  std::sort(evals_gtr.begin(), evals_gtr.end());
  for (size_t i = 0; i < evals_jc.size(); i++) {
    CHECK_LT(fabs(evals_jc[i] - evals_gtr[i]), 0.0001);
  }
  // TODO eigenvectors?
  Eigen::VectorXd frequencies(4);
  frequencies << 0.2, 0.55, 0.1, 0.15;
  Eigen::VectorXd rates(6);
  rates << 1., 0.5, 1., 2., 1., 1.;
  gtr_model->SetParameters({{"frequencies", frequencies}, {"rates", rates}});
  // TODO another eigenvalue/vector test.
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SUBSTITUTION_MODEL_HPP_
