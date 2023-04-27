// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "block_model.hpp"

class SiteModel : public BlockModel {
 public:
  SiteModel(const BlockSpecification::ParamCounts& param_counts)
      : BlockModel(param_counts) {}
  virtual ~SiteModel() = default;

  virtual size_t GetCategoryCount() const = 0;
  virtual const EigenVectorXd& GetCategoryRates() const = 0;
  virtual const EigenVectorXd& GetCategoryProportions() const = 0;
  virtual const EigenVectorXd& GetRateGradient() const = 0;

  static std::unique_ptr<SiteModel> OfSpecification(const std::string& specification);
};

class ConstantSiteModel : public SiteModel {
 public:
  ConstantSiteModel()
      : SiteModel({}), zero_(EigenVectorXd::Zero(1)), one_(EigenVectorXd::Ones(1)) {}

  size_t GetCategoryCount() const override { return 1; }

  const EigenVectorXd& GetCategoryRates() const override { return one_; }

  const EigenVectorXd& GetCategoryProportions() const override { return one_; }

  const EigenVectorXd& GetRateGradient() const override { return zero_; };

  void SetParameters(const EigenVectorXdRef param_vector) override{};

 private:
  EigenVectorXd zero_;
  EigenVectorXd one_;
};

class WeibullSiteModel : public SiteModel {
 public:
  explicit WeibullSiteModel(size_t category_count, double shape)
      : SiteModel({{shape_key_, 1}}),
        category_count_(category_count),
        shape_(shape),
        rate_derivatives_(category_count) {
    category_rates_.resize(category_count);
    category_proportions_.resize(category_count);
    for (size_t i = 0; i < category_count; i++) {
      category_proportions_[i] = 1.0 / category_count;
    }
    UpdateRates();
  }

  size_t GetCategoryCount() const override;
  const EigenVectorXd& GetCategoryRates() const override;
  const EigenVectorXd& GetCategoryProportions() const override;
  const EigenVectorXd& GetRateGradient() const override;

  void SetParameters(const EigenVectorXdRef param_vector) override;

  inline const static std::string shape_key_ = "Weibull_shape";

 private:
  void UpdateRates();

  size_t category_count_;
  double shape_;  // shape of the Weibull distribution
  EigenVectorXd rate_derivatives_;
  EigenVectorXd category_rates_;
  EigenVectorXd category_proportions_;
};
