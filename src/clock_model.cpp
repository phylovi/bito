// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "clock_model.hpp"

std::unique_ptr<ClockModel> ClockModel::OfSpecification(
    const std::string &specification) {
  if (specification == "strict") {
    return std::make_unique<StrictClockModel>();
  }  // else
  Failwith("Clock model not known: " + specification);
}
