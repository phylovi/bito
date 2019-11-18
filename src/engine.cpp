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

std::vector<double> Engine::LogLikelihoods(
    const TreeCollection &tree_collection, EigenMatrixXd &phylo_model_params) {
  return FatBeagleParallelize<double>(FatBeagle::StaticLogLikelihood,
                                      fat_beagles_, tree_collection,
                                      phylo_model_params);
}

std::vector<std::pair<double, std::vector<double>>> Engine::BranchGradients(
    const TreeCollection &tree_collection, EigenMatrixXd &phylo_model_params) {
  return FatBeagleParallelize<std::pair<double, std::vector<double>>>(
      FatBeagle::StaticBranchGradient, fat_beagles_, tree_collection,
      phylo_model_params);
}

