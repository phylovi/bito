// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "site_model.hpp"

#include <cmath>

#include "sugar.hpp"

std::unique_ptr<SiteModel> SiteModel::OfSpecification(
    const std::string& specification) {
  if (specification == "constant") {
    return std::make_unique<ConstantSiteModel>();
  }  // else
  if (specification.rfind("weibull", 0) == 0) {
    size_t category_count = 4;
    auto index = specification.find("+");
    if (index != std::string::npos) {
      std::string cat_count_str = specification.substr(index + 1);
      category_count = stoi(cat_count_str);
    }
    return std::make_unique<WeibullSiteModel>(category_count, 1.0);
  }  // else
  Failwith("Site model not known: " + specification);
}

void WeibullSiteModel::SetParameters(const EigenVectorXdRef param_vector) {
  GetBlockSpecification().CheckParameterVectorSize(param_vector);
  EigenVectorXd shape = ExtractSegment(param_vector, shape_key_);
  shape_ = shape[0];
  UpdateRates();
}

// Discretized Weibull distribution using the median approximation
// Equivalent to the discretized gamma method in Yang 1994.
// The scale (lambda) is fixed to 1
void WeibullSiteModel::UpdateRates() {
  double mean_rate = 0;
  double mean_rate_derivative = 0;
  std::vector<double> deriv_unscaled_rates(category_count_);
  for (size_t i = 0; i < category_count_; i++) {
    double quantile = (2.0 * i + 1.0) / (2.0 * category_count_);
    // Set rate to inverse CDF at quantile.
    category_rates_[i] = pow(-std::log(1.0 - quantile), 1.0 / shape_);
    mean_rate += category_rates_[i];
    // Derivative of unormalized rate i wrt shape
    deriv_unscaled_rates[i] =
        -category_rates_[i] * std::log(-std::log(1.0 - quantile)) / (shape_ * shape_);
    mean_rate_derivative += deriv_unscaled_rates[i];
  }
  mean_rate /= category_count_;
  mean_rate_derivative /= category_count_;

  for (int i = 0; i < category_count_; i++) {
    // Derivative of rate i wrt shape
    // dr_i/dshape = d(ur_i/mean)/dshape = (dur_i* mean - ur_i*dmean)/mean^2
    rate_derivatives_[i] = (deriv_unscaled_rates[i] * mean_rate -
                            category_rates_[i] * mean_rate_derivative) /
                           (mean_rate * mean_rate);
    category_rates_[i] /= mean_rate;
  }
}

size_t WeibullSiteModel::GetCategoryCount() const { return category_count_; }

const EigenVectorXd& WeibullSiteModel::GetCategoryRates() const {
  return category_rates_;
}

const EigenVectorXd& WeibullSiteModel::GetCategoryProportions() const {
  return category_proportions_;
}

const EigenVectorXd& WeibullSiteModel::GetRateGradient() const {
  return rate_derivatives_;
};
