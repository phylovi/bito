// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "engine.hpp"
#include <numeric>
#include "beagle_flag_names.hpp"

Engine::Engine(const EngineSpecification &engine_specification,
               const PhyloModelSpecification &model_specification,
               SitePattern site_pattern)
    : site_pattern_(std::move(site_pattern)) {
  if (engine_specification.thread_count_ == 0) {
    Failwith("Thread count needs to be strictly positive.");
  }  // else
  const auto beagle_preference_flags =
      engine_specification.beagle_flag_vector_.empty()
          ? BEAGLE_FLAG_VECTOR_SSE  // Default flags.
          : std::accumulate(engine_specification.beagle_flag_vector_.begin(),
                            engine_specification.beagle_flag_vector_.end(), 0,
                            std::bit_or<FatBeagle::PackedBeagleFlags>());
  for (size_t i = 0; i < engine_specification.thread_count_; i++) {
    fat_beagles_.push_back(std::make_unique<FatBeagle>(
        model_specification, site_pattern_, beagle_preference_flags,
        engine_specification.use_tip_states_));
  }
  if (!engine_specification.beagle_flag_vector_.empty()) {
    std::cout << "We asked BEAGLE for: "
              << BeagleFlagNames::OfBeagleFlags(beagle_preference_flags) << std::endl;
    auto beagle_flags = fat_beagles_[0]->GetBeagleFlags();
    std::cout << "BEAGLE gave us: " << BeagleFlagNames::OfBeagleFlags(beagle_flags)
              << std::endl;
    if (beagle_flags & BEAGLE_FLAG_PROCESSOR_GPU) {
      std::cout << R"raw(
 ____    ____    __  __      __    __  ______   __    __  __
/\  _`\ /\  _`\ /\ \/\ \    /\ \  /\ \/\  _  \ /\ \  /\ \/\ \
\ \ \L\_\ \ \L\ \ \ \ \ \   \ `\`\\/'/\ \ \L\ \\ `\`\\/'/\ \ \
 \ \ \L_L\ \ ,__/\ \ \ \ \   `\ `\ /'  \ \  __ \`\ `\ /'  \ \ \
  \ \ \/, \ \ \/  \ \ \_\ \    `\ \ \   \ \ \/\ \ `\ \ \   \ \_\
   \ \____/\ \_\   \ \_____\     \ \_\   \ \_\ \_\  \ \_\   \/\_\
    \/___/  \/_/    \/_____/      \/_/    \/_/\/_/   \/_/    \/_/
    )raw";
    }
  }
}

const BlockSpecification &Engine::GetPhyloModelBlockSpecification() const {
  // The BlockSpecification is well defined for an Engine because the interface
  // assures that all of the PhyloModels have the same specification.
  return GetFirstFatBeagle()->GetPhyloModelBlockSpecification();
}

std::vector<double> Engine::LogLikelihoods(const TreeCollection &tree_collection,
                                           const EigenMatrixXdRef phylo_model_params,
                                           const bool rescaling) const {
  return FatBeagleParallelize<double, Tree, TreeCollection>(
      FatBeagle::StaticLogLikelihood, fat_beagles_, tree_collection, phylo_model_params,
      rescaling);
}

std::vector<double> Engine::LogLikelihoods(const RootedTreeCollection &tree_collection,
                                           const EigenMatrixXdRef phylo_model_params,
                                           const bool rescaling) const {
  return FatBeagleParallelize<double, RootedTree, RootedTreeCollection>(
      FatBeagle::StaticRootedLogLikelihood, fat_beagles_, tree_collection,
      phylo_model_params, rescaling);
}

std::vector<std::pair<double, std::vector<double>>> Engine::BranchGradients(
    const TreeCollection &tree_collection, const EigenMatrixXdRef phylo_model_params,
    const bool rescaling) const {
  return FatBeagleParallelize<std::pair<double, std::vector<double>>, Tree,
                              TreeCollection>(FatBeagle::StaticBranchGradient,
                                              fat_beagles_, tree_collection,
                                              phylo_model_params, rescaling);
}

std::vector<std::pair<double, std::vector<double>>> Engine::BranchGradients(
    const RootedTreeCollection &tree_collection,
    const EigenMatrixXdRef phylo_model_params, const bool rescaling) const {
  return FatBeagleParallelize<std::pair<double, std::vector<double>>, RootedTree,
                              RootedTreeCollection>(
      FatBeagle::StaticRootedBranchGradient, fat_beagles_, tree_collection,
      phylo_model_params, rescaling);
}

const FatBeagle *const Engine::GetFirstFatBeagle() const {
  Assert(!fat_beagles_.empty(), "You have no FatBeagles.");
  return fat_beagles_[0].get();
}
