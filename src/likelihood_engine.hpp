// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIKELIHOOD_ENGINE_HPP_
#define SRC_LIKELIHOOD_ENGINE_HPP_

#include <memory>
#include <utility>
#include <vector>
#include "beagle.hpp"
#include "substitution_model.hpp"

class LikelihoodEngine {
 public:
  explicit LikelihoodEngine(SitePattern site_pattern,
                            SubstitutionModelPtr substitution_model,
                            size_t thread_count)
      : site_pattern_(std::move(site_pattern)),
        substitution_model_(std::move(substitution_model)),
        beagle_instances_(thread_count) {
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

  ~LikelihoodEngine() { FinalizeBeagleInstances(); }

  void FinalizeBeagleInstances();

 private:
  SitePattern site_pattern_;
  SubstitutionModelPtr substitution_model_;
  std::vector<beagle::BeagleInstance> beagle_instances_;
};


#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("LikelihoodEngine") {
  auto substitution_model = std::make_unique<GTRModel>();

  // auto engine = LikelihoodEngine(std::move(substitution_model));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_LIKELIHOOD_ENGINE_HPP_
