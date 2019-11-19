// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "phylo_model.hpp"

PhyloModel::PhyloModel(std::unique_ptr<SubstitutionModel> substitution_model,
                       std::unique_ptr<SiteModel> site_model,
                       std::unique_ptr<ClockModel> clock_model)
    : substitution_model_(std::move(substitution_model)),
      site_model_(std::move(site_model)),
      clock_model_(std::move(clock_model)),
      block_specification_({}) {
  size_t next_available_idx = 0;
  // TODO put these strings somewhere.
  substitution_model_->AddToBlockSpecification(
      "substitution total", next_available_idx, block_specification_);
  site_model_->AddToBlockSpecification("site total", next_available_idx,
                                       block_specification_);
  clock_model_->AddToBlockSpecification("clock total", next_available_idx,
                                        block_specification_);
  block_specification_.Insert(BlockSpecification::entire_key_,
                              {0, next_available_idx});
}

std::unique_ptr<PhyloModel> PhyloModel::OfSpecification(
    const PhyloModelSpecification& specification) {
  return std::make_unique<PhyloModel>(
      SubstitutionModel::OfSpecification(specification.substitution_),
      SiteModel::OfSpecification(specification.site_),
      ClockModel::OfSpecification(specification.clock_));
}

void PhyloModel::SetParameters(const EigenVectorXdRef parameterization) {
  substitution_model_->SetParameters(
      ExtractFromParameterization(parameterization, "substitution total"));
  site_model_->SetParameters(
      ExtractFromParameterization(parameterization, "site total"));
  clock_model_->SetParameters(
      ExtractFromParameterization(parameterization, "clock total"));
}

const BlockSpecification& PhyloModel::GetBlockSpecification() const {
  return block_specification_;
}
