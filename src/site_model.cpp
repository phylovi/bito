// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "site_model.hpp"
#include "sugar.hpp"

#include <cmath>

std::unique_ptr<SiteModel> SiteModel::OfSpecification(
    const std::string& specification) {
  if (specification == "constant") {
    return std::make_unique<ConstantSiteModel>();
  } else if (specification.rfind("weibull", 0) == 0) {
    size_t category_count = 4;
    auto index = specification.find("+");
    if (index != std::string::npos) {
      std::string cat_count_str = specification.substr(index + 1);
      category_count = stoi(cat_count_str);
    }
    return std::make_unique<WeibullSiteModel>(category_count);
  }
  Failwith("Site model not known: " + specification);
}

void WeibullSiteModel::SetParameters(const EigenVectorXdRef param_vector) {
  GetBlockSpecification().CheckParameterVectorSize(param_vector);
  Eigen::VectorXd shape = ExtractSegment(param_vector, shape_key_);
  shape_ = shape[0];
  UpdateRates();
}

// Discretized Weibull distribution using the median approximation
// Equivalent to the discretized gamma method in Yang 1994
// The scale (lambda) is fixed to 1
void WeibullSiteModel::UpdateRates() {
  double mean = 0;
  for (size_t i = 0; i < category_count_; i++) {
    double quantile = (2.0 * i + 1.0) / (2.0 * category_count_);
    // Set rate to inverse CDF at quantile
    category_rates_[i] = pow(-std::log(1.0 - quantile), 1.0 / shape_);
    mean += category_rates_[i];
  }
  mean /= category_count_;
  for (int i = 0; i < category_count_; i++) {
    category_rates_[i] /= mean;
  }
}

size_t WeibullSiteModel::GetCategoryCount() { return category_count_; }

const Eigen::VectorXd& WeibullSiteModel::GetCategoryRates() {
  return category_rates_;
}

const Eigen::VectorXd& WeibullSiteModel::GetCategoryProportions() {
  return category_proportions_;
}
