// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Calculation of the ratio and root height gradient, adpated from BEAST.
// https://github.com/beast-dev/beast-mcmc
// Credit to Xiang Ji and Marc Suchard.

#pragma once

#include <numeric>

#include "rooted_tree.hpp"

// \partial{L}/\partial{t_k} = \sum_j \partial{L}/\partial{b_j}
// \partial{b_j}/\partial{t_k}
std::vector<double> HeightGradient(const RootedTree &tree,
                                   const std::vector<double> &branch_gradient);

double GetNodePartial(size_t node_id, size_t leaf_count,
                      const std::vector<double> &heights,
                      const std::vector<double> &ratios,
                      const std::vector<double> &bounds);

// Calculate \partial{t_j}/\partial{r_k}
double GetEpochGradientAddition(
    size_t node_id, size_t child_id, size_t leaf_count,
    const std::vector<double> &heights, const std::vector<double> &ratios,
    const std::vector<double> &bounds,
    const std::vector<double> &ratiosGradientUnweightedLogDensity);

std::vector<double> GetLogTimeArray(const RootedTree &tree);

// Update ratio gradient with \partial{t_j}/\partial{r_k}
std::vector<double> UpdateGradientUnWeightedLogDensity(
    const RootedTree &tree, const std::vector<double> &gradient_height);

double UpdateHeightParameterGradientUnweightedLogDensity(
    const RootedTree &tree, const std::vector<double> &gradient);

std::vector<double> RatioGradientOfHeightGradient(
    const RootedTree &tree, const std::vector<double> &height_gradient);

std::vector<double> RatioGradientOfBranchGradient(
    const RootedTree &tree, const std::vector<double> &branch_gradient);

// This should go away with #205.
EigenVectorXd RatioGradientOfHeightGradientEigen(const RootedTree &tree,
                                                 EigenConstVectorXdRef height_gradient);
