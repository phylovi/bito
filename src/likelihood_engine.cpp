// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "likelihood_engine.hpp"

void LikelihoodEngine::FinalizeBeagleInstances() {
  for (const auto &beagle_instance : beagle_instances_) {
    Assert(beagleFinalizeInstance(beagle_instance) == 0,
           "beagleFinalizeInstance gave nonzero return value!");
  }
  beagle_instances_.clear();
}

