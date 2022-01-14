// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

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
