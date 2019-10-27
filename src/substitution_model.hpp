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
  GTRModel(std::vector<double> rates, std::vector<double> frequencies)
      : rates_(rates) {
    frequencies_ = frequencies;
    need_update_ = true;
    evec_.resize(16);
    ivec_.resize(16);
    eval_.resize(4);
  }
  typedef Eigen::Matrix<double, 4, 4, Eigen::RowMajor> EigenMatrix4d;
  typedef Eigen::Vector4d EigenVector4d;

 protected:
  virtual void UpdateEigenDecomposition();

 private:
  std::vector<double> rates_;

  EigenMatrix4d sqrtPi_;
  EigenMatrix4d sqrtPiInv_;
  EigenMatrix4d qmatrix_;
  EigenMatrix4d eigenvectors_;
  EigenMatrix4d inverse_eigenvectors_;
  EigenVector4d eigenvalues_;
  bool need_update_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("SubstitutionModel") {
  std::vector<double> rates = {1, 1, 1, 1, 1, 1};
  std::vector<double> frequencies = {0.25, 0.25, 0.25, 0.25};
  auto substitution_model = std::make_unique<GTRModel>(rates, frequencies);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SUBSTITUTION_MODEL_HPP_
