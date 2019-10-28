// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SUBSTITUTION_MODEL_HPP_
#define SRC_SUBSTITUTION_MODEL_HPP_
#include <Eigen/Dense>
#include <memory>
#include <vector>

class SubstitutionModel {
 public:
  virtual ~SubstitutionModel() {}

  size_t GetStateCount() const { return frequencies_.size(); }

  const std::vector<double>& GetFrequencies() { return frequencies_; }

  const std::vector<double>& GetEigenVectors() {
    UpdateEigenDecomposition();
    return evec_;
  };
  const std::vector<double>& GetInverseEigenVectors() {
    UpdateEigenDecomposition();
    return ivec_;
  }
  const std::vector<double>& GetEigenValues() {
    UpdateEigenDecomposition();
    return eval_;
  }

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

 protected:
  virtual void UpdateEigenDecomposition() {}
};

class GTRModel : public SubstitutionModel {
 public:
  GTRModel() {
    need_update_ = true;
    evec_.resize(16);
    ivec_.resize(16);
    eval_.resize(4);
    rates_.assign(6, 1.0);
    frequencies_.assign(4, 0.25);
  }

  GTRModel(std::vector<double> rates, std::vector<double> frequencies)
      : GTRModel() {
    rates_ = rates;
    frequencies_ = frequencies;
  }
  typedef Eigen::Matrix<double, 4, 4, Eigen::RowMajor> EigenMatrix4d;
  typedef Eigen::Vector4d EigenVector4d;

 protected:
  virtual void UpdateEigenDecomposition();

 private:
  std::vector<double> rates_;
  bool need_update_;
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
