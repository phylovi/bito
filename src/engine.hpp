// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIKELIHOOD_ENGINE_HPP_
#define SRC_LIKELIHOOD_ENGINE_HPP_

#include <memory>
#include <utility>
#include <vector>
#include "beagle.hpp"
#include "substitution_model.hpp"

class Engine {
 public:
  Engine(SitePattern site_pattern, SubstitutionModelPtr substitution_model,
         size_t thread_count);
  ~Engine();

  std::vector<beagle::BeagleInstance> BeagleInstances() const {
    return beagle_instances_;
  }

  void FinalizeBeagleInstances();

 private:
  SitePattern site_pattern_;
  SubstitutionModelPtr substitution_model_;
  std::vector<beagle::BeagleInstance> beagle_instances_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Engine") {
  auto substitution_model = std::make_unique<GTRModel>();

  // auto engine = Engine(std::move(substitution_model));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_LIKELIHOOD_ENGINE_HPP_
