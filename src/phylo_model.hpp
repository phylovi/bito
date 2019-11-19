// Copyright 2019 libsbn project contributors.
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

class PhyloModel {
 public:
  PhyloModel(std::unique_ptr<SubstitutionModel> substitution_model,
             std::unique_ptr<SiteModel> site_model,
             std::unique_ptr<ClockModel> clock_model);

  SubstitutionModel* GetSubstitutionModel() const {
    return substitution_model_.get();
  }
  SiteModel* GetSiteModel() const { return site_model_.get(); }
  ClockModel* GetClockModel() const { return clock_model_.get(); }

  static std::unique_ptr<PhyloModel> OfSpecification(
      const PhyloModelSpecification& specification);
  // TODO const version?
  void SetParameters(EigenVectorRef parameterization);

  const BlockSpecification& GetBlockSpecification() const;

  // TODO move to cpp file
  // TODO const version? Factor out?
  EigenVectorRef ExtractFromParameterization(EigenVectorRef parameterization,
                                             std::string key) {
    auto [start_idx, parameter_count] = block_specification_.Find(key);
    if (start_idx + parameter_count > parameterization.size()) {
      Failwith("Model parameter " + key +
               " request too long for a parameterization of length " +
               std::to_string(parameterization.size()) + ".");
    }  // else
    return parameterization.segment(start_idx, parameter_count);
  }

 private:
  std::unique_ptr<SubstitutionModel> substitution_model_;
  std::unique_ptr<SiteModel> site_model_;
  std::unique_ptr<ClockModel> clock_model_;

  BlockSpecification block_specification_;
};

#endif  // SRC_PHYLO_MODEL_HPP_
