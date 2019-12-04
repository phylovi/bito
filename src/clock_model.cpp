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

void StrictClockModel::SetParameters(const EigenVectorXdRef param_vector) {
  GetBlockSpecification().CheckParameterVectorSize(param_vector);
  Eigen::VectorXd rate = ExtractSegment(param_vector, rate_key_);
  rate_ = rate[0];
}
