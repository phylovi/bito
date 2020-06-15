// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Calculation of the ratio and root height gradient is adpated from BEAST.
// https://github.com/beast-dev/beast-mcmc
// Credit to Xiang Ji and Marc Suchard.


#include <numeric>
#include "rooted_tree.hpp"

// Calculation of the substitution rate gradient.
// \partial{L}/\partial{r_i} = \partial{L}/\partial{b_i} \partial{b_i}/\partial{r_i}
// For strict clock:
// \partial{L}/\partial{r} = \sum_i \partial{L}/\partial{r_i}
std::vector<double> ClockGradient(const RootedTree &tree,
                                  const std::vector<double> &branch_gradient);

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

std::vector<double> RatioGradientOfBranchGradient(
    const RootedTree &tree, const std::vector<double> &branch_gradient);

std::vector<double> RatioGradientOfHeightGradient(
    const RootedTree &tree, const std::vector<double> &height_gradient);

