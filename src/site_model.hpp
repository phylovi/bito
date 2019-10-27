// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SITE_MODEL_HPP_
#define SRC_SITE_MODEL_HPP_

#include <stdio.h>

#include <vector>

class SiteModelInterface {
 public:
  virtual ~SiteModelInterface() {}

  virtual size_t GetCategoryCount() = 0;

  virtual const std::vector<double>& GetCategoryRates() = 0;

  virtual const std::vector<double>& GetCategoryProportions() = 0;
};

class SiteModel : public SiteModelInterface {
 public:
  SiteModel() : one_(1, 1.0) {}

  virtual size_t GetCategoryCount() { return 1; }

  virtual const std::vector<double>& GetCategoryRates() { return one_; }

  virtual const std::vector<double>& GetCategoryProportions() { return one_; }

 private:
  std::vector<double> one_;
};

class WeibullSiteModel : public SiteModelInterface {
 public:
  WeibullSiteModel(size_t category_count)
      : category_count_(category_count), need_update_(true), shape_(1.0) {
    category_rates_.resize(category_count);
    category_proportions_.assign(category_count, 1.0 / category_count);
  }
  virtual size_t GetCategoryCount() = 0;

  virtual const std::vector<double>& GetCategoryRates() = 0;

  virtual const std::vector<double>& GetCategoryProportions() = 0;

 private:
  void UpdateCategories();

  bool need_update_;
  size_t category_count_;
  double shape_;  // shape of the Weibull distribution
  std::vector<double> category_rates_;
  std::vector<double> category_proportions_;
};
#endif
