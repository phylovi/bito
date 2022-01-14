// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "clock_model.hpp"

std::unique_ptr<ClockModel> ClockModel::OfSpecification(
    const std::string &specification) {
  if (specification == "none") {
    return std::make_unique<NoClockModel>();
  }  // else
  if (specification == "strict") {
    return std::make_unique<StrictClockModel>();
  }  // else
  Failwith("Clock model not known: " + specification);
}

void StrictClockModel::SetParameters(const EigenVectorXdRef param_vector) {
  GetBlockSpecification().CheckParameterVectorSize(param_vector);
  EigenVectorXd rate = ExtractSegment(param_vector, rate_key_);
  rate_ = rate[0];
}
