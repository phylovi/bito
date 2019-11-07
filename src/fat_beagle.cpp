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
  if (finalize_result) {
    std::cout << "beagleFinalizeInstance gave nonzero return value!";
    std::terminate();
  }
}

double FatBeagle::LogLikelihood(
    FatBeagle *fat_beagle,
    // TODO
    // const SubstitutionModel::Parameterization &parameterization,
    const Tree &in_tree) {
  // TODO It's possible that this is just a wrapper for calling a FatBeagle
  // likelihood computation method. Same for BranchGradient.
  Failwith("Not implemented.");
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradient(
    FatBeagle *fat_beagle,
    // TODO
    // const SubstitutionModel::Parameterization &parameterization,
    const Tree &in_tree) {
  Failwith("Not implemented.");
}
