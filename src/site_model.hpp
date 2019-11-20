// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SITE_MODEL_HPP_
#define SRC_SITE_MODEL_HPP_

#include <memory>
#include <vector>
#include "block_model.hpp"

class SiteModel : public BlockModel {
 public:
  SiteModel(const BlockSpecification::ParamCounts& param_counts)
      : BlockModel(param_counts) {}
  virtual ~SiteModel() = default;

  virtual size_t GetCategoryCount() = 0;
  virtual const std::vector<double>& GetCategoryRates() = 0;
  virtual const std::vector<double>& GetCategoryProportions() = 0;

  static std::unique_ptr<SiteModel> OfSpecification(
      const std::string& specification);
};

class ConstantSiteModel : public SiteModel {
 public:
  ConstantSiteModel() : SiteModel({}), one_(1, 1.0) {}

  size_t GetCategoryCount() override { return 1; }

  const std::vector<double>& GetCategoryRates() override { return one_; }

  const std::vector<double>& GetCategoryProportions() override { return one_; }

  void SetParameters(const EigenVectorXdRef parameters) override{};

 private:
  std::vector<double> one_;
};

class WeibullSiteModel : public SiteModel {
 public:
  explicit WeibullSiteModel(size_t category_count)
      : SiteModel({{rates_key_, category_count},
                   {proportions_key_, category_count},
                   {shape_key_, 1}}),
        category_count_(category_count),
        need_update_(true),
        shape_(1.0) {
    category_rates_.resize(category_count);
    category_proportions_.assign(category_count, 1.0 / category_count);
  }

  size_t GetCategoryCount() override;
  const std::vector<double>& GetCategoryRates() override;
  const std::vector<double>& GetCategoryProportions() override;

  void SetParameters(const EigenVectorXdRef parameters) override {
    Failwith("not impelemented");
  };

  inline const static std::string rates_key_ = "Weibull category rates";
  inline const static std::string proportions_key_ =
      "Weibull category proportions";
  inline const static std::string shape_key_ = "Weibull shape";

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
