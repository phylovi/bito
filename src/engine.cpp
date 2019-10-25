// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "engine.hpp"

Engine::Engine(SitePattern site_pattern,
               SubstitutionModelPtr substitution_model, size_t thread_count)
    : site_pattern_(std::move(site_pattern)),
      substitution_model_(std::move(substitution_model)),
      beagle_instances_(thread_count) {
  if (thread_count == 0) {
    Failwith("Thread count needs to be strictly positive.");
  }
  // TODO: allow other subs models
  Assert(substitution_model_->Details().type == SubstitutionModelType::JC,
         "Only JC model allowed.");
  auto make_beagle_instance = [& site_pattern =
                                   std::as_const(this->site_pattern_)]() {
    auto beagle_instance = beagle::CreateInstance(site_pattern);
    beagle::SetJCModel(beagle_instance);
    beagle::PrepareBeagleInstance(beagle_instance, site_pattern);
    return beagle_instance;
  };
  std::generate(beagle_instances_.begin(), beagle_instances_.end(),
                make_beagle_instance);
}

Engine::~Engine() {
  for (const auto& beagle_instance : beagle_instances_) {
    auto finalize_result = beagleFinalizeInstance(beagle_instance);
    if (finalize_result) {
      std::cout << "beagleFinalizeInstance gave nonzero return value!";
      std::terminate();
    }
  }
}

std::vector<double> Engine::LogLikelihoods(
    const TreeCollection& tree_collection, bool rescaling) {
  return beagle::LogLikelihoods(beagle_instances_, tree_collection, rescaling);
}

std::vector<std::pair<double, std::vector<double>>> Engine::BranchGradients(
    const TreeCollection& tree_collection, bool rescaling) {
  return beagle::BranchGradients(beagle_instances_, tree_collection, rescaling);
}

