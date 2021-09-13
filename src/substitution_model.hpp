// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SUBSTITUTION_MODEL_HPP_
#define SRC_SUBSTITUTION_MODEL_HPP_
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "block_model.hpp"
#include "sugar.hpp"

class SubstitutionModel : public BlockModel {
 public:
  SubstitutionModel(const BlockSpecification::ParamCounts& param_counts)
      : BlockModel(param_counts) {}
  virtual ~SubstitutionModel() = default;

  size_t GetStateCount() const { return frequencies_.size(); }

  const EigenMatrixXd& GetQMatrix() const { return Q_; }
  const EigenVectorXd& GetFrequencies() const { return frequencies_; }
  const EigenVectorXd& GetRates() const { return rates_; }
  // We follow BEAGLE in terminology. "Inverse Eigenvectors" means the inverse
  // of the matrix containing the eigenvectors.
  const EigenMatrixXd& GetEigenvectors() const { return eigenvectors_; }
  const EigenMatrixXd& GetInverseEigenvectors() const { return inverse_eigenvectors_; }
  const EigenVectorXd& GetEigenvalues() const { return eigenvalues_; }

  virtual void SetParameters(const EigenVectorXdRef param_vector) = 0;

  // This is the factory method that will be the typical way of buiding
  // substitution models.
  static std::unique_ptr<SubstitutionModel> OfSpecification(
      const std::string& specification);

  inline const static std::string rates_key_ = "substitution model rates";
  inline const static std::string frequencies_key_ = "substitution model frequencies";

 protected:
  EigenVectorXd frequencies_;
  EigenVectorXd rates_;
  EigenMatrixXd eigenvectors_;
  EigenMatrixXd inverse_eigenvectors_;
  EigenVectorXd eigenvalues_;
  EigenMatrixXd Q_;
};

class DNAModel : public SubstitutionModel {
 public:
  DNAModel(const BlockSpecification::ParamCounts& param_counts)
      : SubstitutionModel(param_counts) {
    frequencies_.resize(4);
    eigenvectors_.resize(4, 4);
    inverse_eigenvectors_.resize(4, 4);
    eigenvalues_.resize(4);
    Q_.resize(4, 4);
  }

 protected:
  virtual void UpdateEigendecomposition();
  virtual void UpdateQMatrix() = 0;
  void Update();
};

class JC69Model : public DNAModel {
 public:
  JC69Model() : DNAModel({}) {
    frequencies_ << 0.25, 0.25, 0.25, 0.25;
    Update();
  }

  // No parameters to set for JC!
  void SetParameters(const EigenVectorXdRef) override{};  // NOLINT
  virtual void UpdateEigendecomposition() override;
  void UpdateQMatrix() override;
};

class GTRModel : public DNAModel {
 public:
  explicit GTRModel() : DNAModel({{rates_key_, 6}, {frequencies_key_, 4}}) {
    rates_.resize(6);
    rates_.setConstant(1.0 / 6.0);
    frequencies_ << 0.25, 0.25, 0.25, 0.25;
    Update();
  }
  void SetParameters(const EigenVectorXdRef param_vector) override;

 protected:
  void UpdateQMatrix() override;
};

// The Hasegawa, Kishino and Yano (HKY) susbtitution model.
//
// Reference:
// Hasegawa, M., Kishino, H. and Yano, T.A., 1985. Dating of the human-ape splitting
// by a molecular clock of mitochondrial DNA. Journal of molecular evolution, 22(2),
// pp.160-174.
class HKYModel : public DNAModel {
 public:
  explicit HKYModel() : DNAModel({{rates_key_, 1}, {frequencies_key_, 4}}) {
    rates_.resize(1);
    rates_.setConstant(1.0);
    frequencies_ << 0.25, 0.25, 0.25, 0.25;
    Update();
  }
  void SetParameters(const EigenVectorXdRef param_vector) override;

 protected:
  virtual void UpdateEigendecomposition() override;
  void UpdateQMatrix() override;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <algorithm>
TEST_CASE("SubstitutionModel") {
  auto CheckEigenvalueEquality = [](EigenVectorXd eval1, EigenVectorXd eval2) {
    std::sort(eval1.begin(), eval1.end());
    std::sort(eval2.begin(), eval2.end());
    CheckVectorXdEquality(eval1, eval2, 0.0001);
  };
  auto gtr_model = std::make_unique<GTRModel>();
  auto hky_model = std::make_unique<HKYModel>();
  auto jc_model = std::make_unique<JC69Model>();
  // Test 1: First we test using the "built in" default values.
  CheckEigenvalueEquality(jc_model->GetEigenvalues(), gtr_model->GetEigenvalues());
  CheckEigenvalueEquality(jc_model->GetEigenvalues(), hky_model->GetEigenvalues());
  EigenVectorXd param_vector(10);
  // Test 2: Now try out ParameterSegmentMapOf.
  gtr_model = std::make_unique<GTRModel>();
  // First zero out our param_vector.
  param_vector.setZero();
  // We can use ParameterSegmentMapOf to get two "views" into our parameter
  // vector.
  auto parameter_map =
      gtr_model->GetBlockSpecification().ParameterSegmentMapOf(param_vector);
  auto frequencies = parameter_map.at(SubstitutionModel::frequencies_key_);
  auto rates = parameter_map.at(SubstitutionModel::rates_key_);
  // When we modify the contents of these views, that changes param_vector.
  frequencies.setConstant(0.25);
  rates.setConstant(1.0 / 6.0);
  // We can then set param_vector and go forward as before.
  gtr_model->SetParameters(param_vector);
  CheckEigenvalueEquality(jc_model->GetEigenvalues(), gtr_model->GetEigenvalues());
  // Test 3: Compare to eigenvalues from R.
  frequencies << 0.479367, 0.172572, 0.140933, 0.207128;
  rates << 0.060602, 0.402732, 0.028230, 0.047910, 0.407249, 0.053277;
  gtr_model->SetParameters(param_vector);
  EigenVectorXd eigen_values_r(4);
  eigen_values_r << -2.567992e+00, -1.760838e+00, -4.214918e-01, 1.665335e-16;
  CheckEigenvalueEquality(eigen_values_r, gtr_model->GetEigenvalues());
  // Test HKY against GTR
  EigenVectorXd hky_param_vector(5);
  hky_param_vector.setZero();
  auto hky_parameter_map =
      hky_model->GetBlockSpecification().ParameterSegmentMapOf(hky_param_vector);
  auto hky_frequencies = hky_parameter_map.at(SubstitutionModel::frequencies_key_);
  auto hky_kappa = hky_parameter_map.at(SubstitutionModel::rates_key_);
  hky_frequencies << 0.1, 0.2, 0.3, 0.4;
  hky_kappa.setConstant(3.0);
  hky_model->SetParameters(hky_param_vector);
  frequencies << 0.1, 0.2, 0.3, 0.4;
  rates << 0.1, 0.3, 0.1, 0.1, 0.3, 0.1;
  gtr_model->SetParameters(param_vector);
  CheckEigenvalueEquality(gtr_model->GetEigenvalues(), hky_model->GetEigenvalues());
  CHECK(gtr_model->GetQMatrix().isApprox(hky_model->GetQMatrix()));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SUBSTITUTION_MODEL_HPP_
