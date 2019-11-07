// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_PHYLO_MODEL_HPP_
#define SRC_PHYLO_MODEL_HPP_

#include <memory>

#include "clock_model.hpp"
#include "site_model.hpp"
#include "substitution_model.hpp"

class PhyloModel {
 public:
  PhyloModel(std::unique_ptr<SubstitutionModel> substitution_model,
             std::unique_ptr<SiteModel> site_model,
             std::unique_ptr<ClockModel> clock_model);

 private:
  std::unique_ptr<SubstitutionModel> substitution_model_;
  std::unique_ptr<SiteModel> site_model_;
  std::unique_ptr<ClockModel> clock_model_;
};


#endif  // SRC_PHYLO_MODEL_HPP_
