// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_GRADIENT_HPP_
#define SRC_TREE_GRADIENT_HPP_

#include <map>
#include <vector>

using GradientMap = std::map<std::string, std::vector<double>>;

struct PhyloGradient {
  PhyloGradient() = default;
  PhyloGradient(double log_likelihood, GradientMap& gradient)
      : log_likelihood_(log_likelihood), gradient_(gradient){};

  double log_likelihood_;
  GradientMap gradient_;
};

#endif  // SRC_TREE_GRADIENT_HPP_
