// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ENGINE_HPP_
#define SRC_ENGINE_HPP_

#include <memory>
#include <utility>
#include <vector>
#include "beagle.hpp"
#include "substitution_model.hpp"

// "Engine" is short for "phylogenetic likelihood computation engine".
class Engine {
 public:
  Engine(SitePattern site_pattern, SubstitutionModelPtr substitution_model,
         size_t thread_count);
  ~Engine();
  // Delete (copy + move) x (constructor + assignment)
  Engine(const Engine &) = delete;
  Engine(const Engine &&) = delete;
  Engine &operator=(const Engine &) = delete;
  Engine &operator=(const Engine &&) = delete;

  std::vector<beagle::BeagleInstance> BeagleInstances() const {
    return beagle_instances_;
  }

  std::vector<double> LogLikelihoods(const TreeCollection &tree_collection,
                                     bool rescaling);
  std::vector<std::pair<double, std::vector<double>>> BranchGradients(
      const TreeCollection &tree_collection, bool rescaling);

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

#endif  // SRC_ENGINE_HPP_
