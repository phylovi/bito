// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "engine.hpp"

Engine::Engine(const PhyloModelSpecification &specification,
               SitePattern site_pattern, size_t thread_count)
    : site_pattern_(std::move(site_pattern)) {
  if (thread_count == 0) {
    Failwith("Thread count needs to be strictly positive.");
  }
  for (size_t i = 0; i < thread_count; i++) {
    fat_beagles_.push_back(
        std::make_unique<FatBeagle>(specification, site_pattern_));
  }
}

const BlockSpecification &Engine::GetPhyloModelBlockSpecification() const {
  // The BlockSpecification is well defined for an Engine because the interface
  // assures that all of the PhyloModels have the same specification.
  return GetFirstFatBeagle()->GetPhyloModelBlockSpecification();
}

std::vector<double> Engine::LogLikelihoods(
    const TreeCollection &tree_collection,
    const EigenMatrixXdRef phylo_model_params) const {
  return FatBeagleParallelize<double>(FatBeagle::StaticLogLikelihood,
                                      fat_beagles_, tree_collection,
                                      phylo_model_params);
}

std::vector<std::pair<double, std::vector<double>>> Engine::BranchGradients(
    const TreeCollection &tree_collection,
    const EigenMatrixXdRef phylo_model_params) const {
  return FatBeagleParallelize<std::pair<double, std::vector<double>>>(
      FatBeagle::StaticBranchGradient, fat_beagles_, tree_collection,
      phylo_model_params);
}

const FatBeagle *const Engine::GetFirstFatBeagle() const {
  Assert(!fat_beagles_.empty(), "You have no FatBeagles.");
  return fat_beagles_[0].get();
}
