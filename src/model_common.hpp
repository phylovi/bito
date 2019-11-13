// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_MODEL_COMMON_HPP_
#define SRC_MODEL_COMMON_HPP_

#include <Eigen/Dense>
#include <unordered_map>

using ModelParameterization = std::unordered_map<std::string, Eigen::VectorXd>;

#endif  // SRC_MODEL_COMMON_HPP_
