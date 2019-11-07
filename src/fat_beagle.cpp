// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "fat_beagle.hpp"
#include "beagle.hpp"

FatBeagle::FatBeagle(const PhyloModel &phylo_model,
                     const SitePattern &site_pattern)
    : phylo_model_(phylo_model), rescaling_(false) {
  beagle_instance_ = beagle::CreateInstance(site_pattern);
  beagle::PrepareBeagleInstance(beagle_instance_, site_pattern);
  beagle::SetJCModel(beagle_instance_);
};

FatBeagle::~FatBeagle() {
  auto finalize_result = beagleFinalizeInstance(beagle_instance_);
  if (finalize_result != 0) {
    std::cout << "beagleFinalizeInstance gave nonzero return value!";
    std::terminate();
  }
}

double FatBeagle::LogLikelihood(
    const Tree &in_tree) {
  return beagle::LogLikelihood(beagle_instance_, in_tree, rescaling_);
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradient(
    const Tree &in_tree) {
  return beagle::BranchGradient(beagle_instance_, in_tree, rescaling_);
}

double FatBeagle::StaticLogLikelihood(
    FatBeagle *fat_beagle,
    const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "Null FatBeagle pointer!");
  return fat_beagle->LogLikelihood(in_tree);
}

std::pair<double, std::vector<double>> FatBeagle::StaticBranchGradient(
    FatBeagle *fat_beagle,
    const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "Null FatBeagle pointer!");
  return fat_beagle->BranchGradient(in_tree);
}
