// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "fat_beagle.hpp"
#include "beagle.hpp"

FatBeagle::FatBeagle(const PhyloModel &phylo_model,
                     const SitePattern &site_pattern)
    : phylo_model_(phylo_model), rescaling_(false) {
  beagle_instance_ = beagle::CreateInstance(site_pattern);
  beagle::PrepareBeagleInstance(beagle_instance_, site_pattern);
};

FatBeagle::~FatBeagle() {
  auto finalize_result = beagleFinalizeInstance(beagle_instance_);
  if (finalize_result != 0) {
    std::cout << "beagleFinalizeInstance gave nonzero return value!";
    std::terminate();
  }
}

double FatBeagle::LogLikelihood(
    // TODO
    // const SubstitutionModel::Parameterization &parameterization,
    const Tree &in_tree) {
  // TODO It's possible that this is just a wrapper for calling a FatBeagle
  // likelihood computation method. Same for BranchGradient.
  Failwith("Not implemented.");
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradient(
    // TODO
    // const SubstitutionModel::Parameterization &parameterization,
    const Tree &in_tree) {
  Failwith("Not implemented.");
}

double FatBeagle::StaticLogLikelihood(
    FatBeagle *fat_beagle,
    // TODO
    // const SubstitutionModel::Parameterization &parameterization,
    const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "Null FatBeagle pointer!");
  return fat_beagle->LogLikelihood(in_tree);
}

std::pair<double, std::vector<double>> FatBeagle::StaticBranchGradient(
    FatBeagle *fat_beagle,
    // TODO
    // const SubstitutionModel::Parameterization &parameterization,
    const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "Null FatBeagle pointer!");
  return fat_beagle->BranchGradient(in_tree);
}
