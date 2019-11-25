// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// "Engine" is short for "phylogenetic likelihood computation engine".
// This engine has FatBeagles as cylinders.

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

  const BlockSpecification &GetPhyloModelBlockSpecification() const;

  std::vector<double> LogLikelihoods(
      const TreeCollection &tree_collection,
      const EigenMatrixXdRef phylo_model_params) const;
  std::vector<std::pair<double, std::vector<double>>> BranchGradients(
      const TreeCollection &tree_collection,
      const EigenMatrixXdRef phylo_model_params) const;

 private:
  SitePattern site_pattern_;
  std::vector<std::unique_ptr<FatBeagle>> fat_beagles_;

  const FatBeagle *const GetFirstFatBeagle() const;
};

#endif  // SRC_ENGINE_HPP_
