// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_GRADIENT_HPP_
#define SRC_TREE_GRADIENT_HPP_

#include <vector>
#include "eigen_sugar.hpp"

struct TreeGradient {
  TreeGradient() = default;
  TreeGradient(double log_likelihood, std::vector<double> site_model_gradient,
               std::vector<double> substitution_model_gradient)
      : log_likelihood_(log_likelihood),
        site_model_(site_model_gradient),
        substitution_model_(substitution_model_gradient){};

  double log_likelihood_;
  std::vector<double> site_model_;
  std::vector<double> substitution_model_;
};

struct UnrootedTreeGradient : TreeGradient {
  UnrootedTreeGradient() = default;
  UnrootedTreeGradient(double log_likelihood,
                       std::vector<double> branch_length_gradient,
                       std::vector<double> site_model_gradient,
                       std::vector<double> substitution_model_gradient)
      : TreeGradient(log_likelihood, site_model_gradient, substitution_model_gradient),
        branch_lengths_(branch_length_gradient){};

  std::vector<double> branch_lengths_;
};

struct RootedTreeGradient : TreeGradient {
  RootedTreeGradient() = default;
  RootedTreeGradient(double log_likelihood, std::vector<double> branch_length_gradient,
                     EigenVectorXd ratios_root_height_gradient,
                     EigenVectorXd clock_model_gradient,
                     EigenVectorXd site_model_gradient,
                     EigenVectorXd substitution_model_gradient)
      : TreeGradient(log_likelihood, site_model_gradient, substitution_model_gradient),
        branch_lengths_(branch_length_gradient),
        clock_model_(clock_model_gradient),
        ratios_root_height_(ratios_root_height_gradient){};

  EigenVectorXd branch_lengths_;
  EigenVectorXd clock_model_;
  EigenVectorXd ratios_root_height_;
};

#endif  // SRC_TREE_GRADIENT_HPP_
