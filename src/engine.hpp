// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ENGINE_HPP_
#define SRC_ENGINE_HPP_

#include <memory>
#include <utility>
#include <vector>
#include "fat_beagle.hpp"
#include "phylo_model.hpp"
#include "site_pattern.hpp"
#include "tree_collection.hpp"

// "Engine" is short for "phylogenetic likelihood computation engine".
class Engine {
 public:
  Engine(PhyloModel phylo_model, SitePattern site_pattern, size_t thread_count);
  ~Engine() = default;
  // Delete (copy + move) x (constructor + assignment)
  Engine(const Engine &) = delete;
  Engine(const Engine &&) = delete;
  Engine &operator=(const Engine &) = delete;
  Engine &operator=(const Engine &&) = delete;


  std::vector<double> LogLikelihoods(const TreeCollection &tree_collection);
  std::vector<std::pair<double, std::vector<double>>> BranchGradients(
      const TreeCollection &tree_collection);

 private:
  PhyloModel phylo_model_;
  SitePattern site_pattern_;
  std::vector<std::unique_ptr<FatBeagle>> fat_beagles_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Engine") {
  auto substitution_model = std::make_unique<GTRModel>();

  // auto engine = Engine(std::move(substitution_model));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ENGINE_HPP_
