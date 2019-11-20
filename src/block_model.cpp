// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "block_model.hpp"
#include "sugar.hpp"

// TODO "parameterization" vs "parameters"?
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

// In the interest of having a single code path, this implementation does twice
// as many lookups as it needs by calling ExtractSegment. My intent is not to
// have it be called frequently-- just once when we're setting things up. So if
// this changes let's do something about it.
BlockModel::ParameterMap BlockModel::ParameterMapOf(
    EigenVectorXdRef parameterization) {
  CheckParametersSize(parameterization);
  ParameterMap parameter_map;
  for (const auto [key, _] : GetBlockSpecification().GetMap()) {
    SafeInsert(parameter_map, key, ExtractSegment(parameterization, key));
  }
  return parameter_map;
}
