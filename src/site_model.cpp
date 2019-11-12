// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "site_model.hpp"
#include "sugar.hpp"

#include <cmath>

std::unique_ptr<SiteModel> SiteModel::OfSpecification(
    const std::string& specification) {
  if (specification == "constant") {
    return std::make_unique<ConstantSiteModel>();
  }  // else
  // TODO Weibull
  Failwith("Site model not known: " + specification);
}

// TODO can you point me to a reference for these computations?
void WeibullSiteModel::UpdateCategories() {
  double mean = 0;
  for (size_t i = 0; i < category_count_; i++) {
    double quantile = (2.0 * i + 1.0) / (2.0 * category_count_);
    // Fix lambda:=1
    category_rates_[i] = pow(-std::log(1.0 - quantile), 1.0 / shape_);
    mean += category_rates_[i];
  }
  mean /= category_count_;
  for (double& rate : category_rates_) {
    rate /= mean;
  }
  need_update_ = false;
}

size_t WeibullSiteModel::GetCategoryCount() { return category_count_; }

const std::vector<double>& WeibullSiteModel::GetCategoryRates() {
  if (need_update_) {
    UpdateCategories();
  }
  return category_rates_;
}

const std::vector<double>& WeibullSiteModel::GetCategoryProportions() {
  if (need_update_) {
    UpdateCategories();
  }
  return category_proportions_;
}
