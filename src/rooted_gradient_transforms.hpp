// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Calculation of the ratio and root height gradient, adpated from BEAST.
// https://github.com/beast-dev/beast-mcmc
// Credit to Xiang Ji and Marc Suchard.

#include <numeric>
#include "rooted_tree.hpp"

// \partial{L}/\partial{t_k} = \sum_j \partial{L}/\partial{b_j}
// \partial{b_j}/\partial{t_k}
EigenVectorXd HeightGradient(const RootedTree &tree,
                             EigenConstVectorXdRef branch_gradient);

double GetNodePartial(size_t node_id, size_t leaf_count, EigenConstVectorXdRef heights,
                      EigenConstVectorXdRef ratios, EigenConstVectorXdRef bounds);

// Calculate \partial{t_j}/\partial{r_k}
double GetEpochGradientAddition(
    size_t node_id, size_t child_id, size_t leaf_count, EigenConstVectorXdRef heights,
    EigenConstVectorXdRef ratios, EigenConstVectorXdRef bounds,
    EigenConstVectorXdRef ratiosGradientUnweightedLogDensity);

EigenVectorXd GetLogTimeArray(const RootedTree &tree);

// Update ratio gradient with \partial{t_j}/\partial{r_k}
EigenVectorXd UpdateGradientUnWeightedLogDensity(const RootedTree &tree,
                                                 EigenConstVectorXdRef gradient_height);

double UpdateHeightParameterGradientUnweightedLogDensity(
    const RootedTree &tree, EigenConstVectorXdRef gradient);

EigenVectorXd RatioGradientOfHeightGradient(const RootedTree &tree,
                                            EigenConstVectorXdRef height_gradient);

EigenVectorXd RatioGradientOfBranchGradient(const RootedTree &tree,
                                            std::vector<double> &branch_gradient);
