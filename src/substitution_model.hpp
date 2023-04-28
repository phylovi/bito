// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "block_model.hpp"
#include "eigen_sugar.hpp"
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

  inline const static std::string rates_key_ = "substitution_model_rates";
  inline const static std::string frequencies_key_ = "substitution_model_frequencies";

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
