// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "phylo_model.hpp"

PhyloModel::PhyloModel(std::unique_ptr<SubstitutionModel> substitution_model,
                       std::unique_ptr<SiteModel> site_model,
                       std::unique_ptr<ClockModel> clock_model)
    : substitution_model_(std::move(substitution_model)),
      site_model_(std::move(site_model)),
      clock_model_(std::move(clock_model)) {}

std::unique_ptr<PhyloModel> PhyloModel::OfSpecification(
    const PhyloModelSpecification& specification) {
  return std::make_unique<PhyloModel>(
      SubstitutionModel::OfSpecification(specification.substitution_),
      SiteModel::OfSpecification(specification.site_),
      ClockModel::OfSpecification(specification.clock_));
}

void PhyloModel::SetParameters(const Eigen::VectorXd& parameterization) {
  substitution_model_->SetParameters(parameterization);
  site_model_->SetParameters(parameterization);
  clock_model_->SetParameters(parameterization);
}

BlockSpecification PhyloModel::GetBlockSpecification() const {
  size_t next_available_idx = 0;
  BlockSpecification blocks;
  substitution_model_->AddToBlockSpecification(next_available_idx, blocks);
  site_model_->AddToBlockSpecification(next_available_idx, blocks);
  clock_model_->AddToBlockSpecification(next_available_idx, blocks);
  return blocks;
}
