// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// "Engine" is short for "phylogenetic likelihood computation engine".

#ifndef SRC_ENGINE_HPP_
#define SRC_ENGINE_HPP_

#include <memory>
#include <utility>
#include <vector>
#include "fat_beagle.hpp"
#include "phylo_model.hpp"
#include "site_pattern.hpp"
#include "tree_collection.hpp"

class Engine {
 public:
  Engine(const PhyloModelSpecification &specification, SitePattern site_pattern,
         size_t thread_count);

  // TODO const
  BlockSpecification GetBlockSpecification() const;

  std::vector<double> LogLikelihoods(const TreeCollection &tree_collection,
                                     const EigenMatrixXdRef phylo_model_params);
  std::vector<std::pair<double, std::vector<double>>> BranchGradients(
      const TreeCollection &tree_collection,
      const EigenMatrixXdRef phylo_model_params);

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
