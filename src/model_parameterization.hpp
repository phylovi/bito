// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_MODEL_PARAMETERIZATION_HPP_
#define SRC_MODEL_PARAMETERIZATION_HPP_

#include <Eigen/Dense>
#include <unordered_map>
#include "sugar.hpp"

using ModelParameterization = std::unordered_map<std::string, Eigen::VectorXd>;
using ModelParameterizationVector = std::vector<ModelParameterization>;

void SetFromParameterization(ModelParameterization parameterization,
                             std::string key,
                             Eigen::VectorXd &parameter_to_set);

#endif  // SRC_MODEL_PARAMETERIZATION_HPP_
