// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_GRADIENT_HPP_
#define SRC_TREE_GRADIENT_HPP_

#include <vector>

struct PhyloGradient {
  PhyloGradient() = default;
  PhyloGradient(double log_likelihood, std::vector<double> site_model_gradient,
                std::vector<double> substitution_model_gradient)
      : log_likelihood_(log_likelihood),
        site_model_(site_model_gradient),
        substitution_model_(substitution_model_gradient){};

  double log_likelihood_;
  std::vector<double> site_model_;
  std::vector<double> substitution_model_;
};

struct UnrootedPhyloGradient : PhyloGradient {
  UnrootedPhyloGradient() = default;
  UnrootedPhyloGradient(double log_likelihood,
                        std::vector<double> branch_length_gradient,
                        std::vector<double> site_model_gradient,
                        std::vector<double> substitution_model_gradient)
      : PhyloGradient(log_likelihood, site_model_gradient, substitution_model_gradient),
        branch_lengths_(branch_length_gradient){};

  std::vector<double> branch_lengths_;
};

struct RootedPhyloGradient : PhyloGradient {
  RootedPhyloGradient() = default;
  RootedPhyloGradient(double log_likelihood, std::vector<double> branch_length_gradient,
                      std::vector<double> ratios_root_height_gradient,
                      std::vector<double> clock_model_gradient,
                      std::vector<double> site_model_gradient,
                      std::vector<double> substitution_model_gradient)
      : PhyloGradient(log_likelihood, site_model_gradient, substitution_model_gradient),
        branch_lengths_(branch_length_gradient),
        clock_model_(clock_model_gradient),
        ratios_root_height_(ratios_root_height_gradient){};

  std::vector<double> branch_lengths_;
  std::vector<double> clock_model_;
  std::vector<double> ratios_root_height_;
};

#endif  // SRC_TREE_GRADIENT_HPP_
