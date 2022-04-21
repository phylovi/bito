// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <memory>

#include "clock_model.hpp"
#include "site_model.hpp"
#include "substitution_model.hpp"
#include "phylo_flags.hpp"

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

  inline const static std::string entire_substitution_key_ = "entire_substitution";
  inline const static std::string entire_site_key_ = "entire_site";
  inline const static std::string entire_clock_key_ = "entire_clock";

 private:
  std::unique_ptr<SubstitutionModel> substitution_model_;
  std::unique_ptr<SiteModel> site_model_;
  std::unique_ptr<ClockModel> clock_model_;
};

// Mapkeys for PhyloModel Parameters
namespace PhyloModelMapkeys {
// Map keys
inline static const auto substitution_model_ =
    PhyloMapkey("SUBSTITUTION_MODEL", PhyloModel::entire_substitution_key_);
inline static const auto substitution_model_rates_ =
    PhyloMapkey("SUBSTITUTION_MODEL_RATES", SubstitutionModel::rates_key_);
inline static const auto substitution_model_frequencies_ =
    PhyloMapkey("SUBSTITUTION_MODEL_FREQUENCIES", SubstitutionModel::frequencies_key_);
inline static const auto site_model =
    PhyloMapkey("SITE_MODEL", PhyloModel::entire_site_key_);
inline static const auto clock_model_ =
    PhyloMapkey("CLOCK_MODEL", PhyloModel::entire_clock_key_);
inline static const auto clock_model_rates_ =
    PhyloMapkey("CLOCK_MODEL_RATES", StrictClockModel::rate_key_);

inline static const auto set_ =
    PhyloMapkeySet("PhyloModel", {substitution_model_, substitution_model_rates_,
                                  substitution_model_frequencies_, site_model,
                                  clock_model_, clock_model_rates_});
}  // namespace PhyloModelMapkeys
