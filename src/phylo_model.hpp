// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_PHYLO_MODEL_HPP_
#define SRC_PHYLO_MODEL_HPP_

#include <memory>

#include "clock_model.hpp"
#include "site_model.hpp"
#include "substitution_model.hpp"

struct PhyloModelSpecification {
  std::string substitution_;
  std::string site_;
  std::string clock_;
};

class PhyloModel : public BlockModel {
 public:
  PhyloModel(std::unique_ptr<SubstitutionModel> substitution_model,
             std::unique_ptr<SiteModel> site_model,
             std::unique_ptr<ClockModel> clock_model);

  SubstitutionModel* GetSubstitutionModel() const { return substitution_model_.get(); }
  SiteModel* GetSiteModel() const { return site_model_.get(); }
  ClockModel* GetClockModel() const { return clock_model_.get(); }

  static std::unique_ptr<PhyloModel> OfSpecification(
      const PhyloModelSpecification& specification);
  void SetParameters(const EigenVectorXdRef param_vector) override;

  inline const static std::string entire_substitution_key_ = "entire substitution";
  inline const static std::string entire_site_key_ = "entire site";
  inline const static std::string entire_clock_key_ = "entire clock";

 private:
  std::unique_ptr<SubstitutionModel> substitution_model_;
  std::unique_ptr<SiteModel> site_model_;
  std::unique_ptr<ClockModel> clock_model_;
};

#endif  // SRC_PHYLO_MODEL_HPP_
