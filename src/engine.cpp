// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "engine.hpp"

Engine::Engine(PhyloModel phylo_model, SitePattern site_pattern,
               size_t thread_count)
    : phylo_model_(std::move(phylo_model)),
      site_pattern_(std::move(site_pattern)) {
  // TODO actually make fat beagle vector
  if (thread_count == 0) {
    Failwith("Thread count needs to be strictly positive.");
  }
  auto make_fat_beagle = [&phylo_model = std::as_const(this->phylo_model_),
                          &site_pattern =
                              std::as_const(this->site_pattern_)]() {
    return std::make_unique<FatBeagle>(phylo_model, site_pattern);
  };
  std::generate(fat_beagles_.begin(), fat_beagles_.end(), make_fat_beagle);
}

std::vector<double> Engine::LogLikelihoods(
    const TreeCollection &tree_collection) {
  return Parallelize<double>(FatBeagle::LogLikelihood, fat_beagles_,
                             tree_collection);
}

std::vector<std::pair<double, std::vector<double>>> Engine::BranchGradients(
    const TreeCollection &tree_collection) {
  return Parallelize<std::pair<double, std::vector<double>>>(
      FatBeagle::BranchGradient, fat_beagles_, tree_collection);
}

