// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "block_model.hpp"
#include "sugar.hpp"

EigenVectorXdRef BlockModel::ExtractSegment(EigenVectorXdRef parameterization,
                                            std::string key) {
  auto [start_idx, parameter_count] = block_specification_.Find(key);
  if (start_idx + parameter_count > parameterization.size()) {
    Failwith("Model parameter " + key +
             " request too long for a parameterization of length " +
             std::to_string(parameterization.size()) + ".");
  }  // else
  return parameterization.segment(start_idx, parameter_count);
}

