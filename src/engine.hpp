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
  Engine(const PhyloModelSpecification &specification, SitePattern site_pattern,
         size_t thread_count);
  ~Engine() = default;
  // Delete (copy + move) x (constructor + assignment)
  Engine(const Engine &) = delete;
  Engine(const Engine &&) = delete;
  Engine &operator=(const Engine &) = delete;
  Engine &operator=(const Engine &&) = delete;

  FatBeagle *GetFatBeagle(size_t idx) {
    Assert(idx < fat_beagles_.size(), "FatBeagle index out of range.");
    return fat_beagles_[idx].get();
  }

  BlockSpecification GetBlockSpecification() const {
    return GetFirstFatBeagle()->GetBlockSpecification();
  }

  std::vector<double> LogLikelihoods(const TreeCollection &tree_collection,
                                     const EigenMatrixXd &phylo_model_params);
  std::vector<std::pair<double, std::vector<double>>> BranchGradients(
      const TreeCollection &tree_collection,
      const EigenMatrixXd &phylo_model_params);

 private:
  SitePattern site_pattern_;
  std::vector<std::unique_ptr<FatBeagle>> fat_beagles_;

  FatBeagle *GetFirstFatBeagle() const {
    Assert(!fat_beagles_.empty(), "You have no FatBeagles.");
    return fat_beagles_[0].get();
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Engine") {
  auto substitution_model = std::make_unique<GTRModel>();

  // auto engine = Engine(std::move(substitution_model));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ENGINE_HPP_
