// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "engine.hpp"

Engine::Engine(PhyloModel phylo_model, SitePattern site_pattern,
               size_t thread_count)
    : phylo_model_(std::move(phylo_model)),
      site_pattern_(std::move(site_pattern)) {
  if (thread_count == 0) {
    Failwith("Thread count needs to be strictly positive.");
  }
  for (size_t i = 0; i < thread_count; i++) {
    fat_beagles_.push_back(
        std::make_unique<FatBeagle>(phylo_model_, site_pattern_));
  }
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

