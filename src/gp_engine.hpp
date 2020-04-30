// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A crappy version of Engine that only does JC, but does do our GP calculations.

#ifndef SRC_GP_ENGINE_HPP_
#define SRC_GP_ENGINE_HPP_

#include "gp_operation.hpp"
#include "site_pattern.hpp"

class GPEngine {
 public:
  GPEngine(SitePattern site_pattern);

 private:
  SitePattern site_pattern_;
};

#endif  // SRC_GP_ENGINE_HPP_
