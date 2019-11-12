// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SITE_MODEL_HPP_
#define SRC_SITE_MODEL_HPP_

#include <stdio.h>

#include <memory>
#include <vector>

class SiteModel {
 public:
  virtual ~SiteModel() = default;

  virtual size_t GetCategoryCount() = 0;
  virtual const std::vector<double>& GetCategoryRates() = 0;
  virtual const std::vector<double>& GetCategoryProportions() = 0;

  static std::unique_ptr<SiteModel> OfSpecification(
      const std::string& specification);
};

class ConstantSiteModel : public SiteModel {
 public:
  ConstantSiteModel() : one_(1, 1.0) {}

  size_t GetCategoryCount() override { return 1; }

  const std::vector<double>& GetCategoryRates() override { return one_; }

  const std::vector<double>& GetCategoryProportions() override { return one_; }

 private:
  std::vector<double> one_;
};

class WeibullSiteModel : public SiteModel {
 public:
  explicit WeibullSiteModel(size_t category_count)
      : category_count_(category_count), need_update_(true), shape_(1.0) {
    category_rates_.resize(category_count);
    category_proportions_.assign(category_count, 1.0 / category_count);
  }
  virtual size_t GetCategoryCount() = 0;

  virtual const std::vector<double>& GetCategoryRates() = 0;

  virtual const std::vector<double>& GetCategoryProportions() = 0;

 private:
  void UpdateCategories();

  size_t category_count_;
  // TODO Would it be possible to recalculate things using setters rather than
  // having an update flag?
  bool need_update_;
  double shape_;  // shape of the Weibull distribution
  std::vector<double> category_rates_;
  std::vector<double> category_proportions_;
};
#endif
