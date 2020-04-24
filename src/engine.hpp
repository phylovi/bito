// Copyright 2019-2020 libsbn project contributors.
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

struct EngineSpecification {
  const size_t thread_count_;
  const std::vector<BeagleFlags> &beagle_flag_vector_;
  const bool use_tip_states_;
};

class Engine {
 public:
  Engine(const EngineSpecification &engine_specification,
         const PhyloModelSpecification &specification, SitePattern site_pattern);

  const BlockSpecification &GetPhyloModelBlockSpecification() const;

  std::vector<double> LogLikelihoods(const TreeCollection &tree_collection,
                                     const EigenMatrixXdRef phylo_model_params,
                                     const bool rescaling) const;
  std::vector<double> LogLikelihoods(const RootedTreeCollection &tree_collection,
                                     const EigenMatrixXdRef phylo_model_params,
                                     const bool rescaling) const;
  std::vector<std::pair<double, std::vector<double>>> BranchGradients(
      const TreeCollection &tree_collection, const EigenMatrixXdRef phylo_model_params,
      const bool rescaling) const;
  std::vector<std::pair<double, std::vector<double>>> BranchGradients(
      const RootedTreeCollection &tree_collection,
      const EigenMatrixXdRef phylo_model_params, const bool rescaling) const;

 private:
  SitePattern site_pattern_;
  std::vector<std::unique_ptr<FatBeagle>> fat_beagles_;

  const FatBeagle *const GetFirstFatBeagle() const;
};

#endif  // SRC_ENGINE_HPP_
