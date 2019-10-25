// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIKELIHOOD_ENGINE_HPP_
#define SRC_LIKELIHOOD_ENGINE_HPP_

#include <memory>
#include <utility>
#include "substitution_model.hpp"

class LikelihoodEngine {
 public:
  explicit LikelihoodEngine(SubstitutionModelPtr substitution_model)
      : substitution_model_(std::move(substitution_model)) {}

 private:
  SubstitutionModelPtr substitution_model_;
};


#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("LikelihoodEngine") {
  auto substitution_model = std::make_unique<GTRModel>();

  auto engine = LikelihoodEngine(std::move(substitution_model));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_LIKELIHOOD_ENGINE_HPP_
