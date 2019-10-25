// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "likelihood_engine.hpp"

LikelihoodEngine::LikelihoodEngine(SitePattern site_pattern,
                                   SubstitutionModelPtr substitution_model,
                                   size_t thread_count)
    : site_pattern_(std::move(site_pattern)),
      substitution_model_(std::move(substitution_model)),
      beagle_instances_(thread_count) {
  // TODO: something else
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

LikelihoodEngine::~LikelihoodEngine() { FinalizeBeagleInstances(); }

void LikelihoodEngine::FinalizeBeagleInstances() {
  for (const auto &beagle_instance : beagle_instances_) {
    Assert(beagleFinalizeInstance(beagle_instance) == 0,
           "beagleFinalizeInstance gave nonzero return value!");
  }
  beagle_instances_.clear();
}

