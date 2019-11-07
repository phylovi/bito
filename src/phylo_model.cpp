// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "phylo_model.hpp"

PhyloModel::PhyloModel(std::unique_ptr<SubstitutionModel> substitution_model,
                       std::unique_ptr<SiteModel> site_model,
                       std::unique_ptr<ClockModel> clock_model)
    : substitution_model_(std::move(substitution_model)),
      site_model_(std::move(site_model)),
      clock_model_(std::move(clock_model)) {}
