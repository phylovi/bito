// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "model_parameterization.hpp"

void SetFromParameterization(ModelParameterization parameterization,
                             std::string key,
                             Eigen::VectorXd &parameter_to_set) {
  auto search = parameterization.find(key);
  if (search == parameterization.end()) {
    // TODO this function does nothing if key is not found in parameterization.
    // This seems error-prone?
    return;
  }  // else
  auto parameter = search->second;
  if (parameter.size() != parameter_to_set.size()) {
    Failwith("Model parameter " + key + " has length " +
             std::to_string(parameter.size()) + " but expected size was " +
             std::to_string(parameter_to_set.size()) + ".");
  }  // else
  parameter_to_set = parameter;
}
