// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Calculation of the ratio and root height gradient, adpated from BEAST.
// https://github.com/beast-dev/beast-mcmc
// Credit to Xiang Ji and Marc Suchard.
//
// Because this is code adapted from elsewhere, at least for the time being the naming
// conventions are a little different: pascalCase is allowed for variables.

#include <numeric>

#include "rooted_tree.hpp"

// \partial{L}/\partial{t_k} = \sum_j \partial{L}/\partial{b_j}
// \partial{b_j}/\partial{t_k}
std::vector<double> HeightGradient(const RootedTree &tree,
                                   const std::vector<double> &branch_gradient) {
  auto root_id = tree.Topology()->Id();
  std::vector<double> height_gradient(tree.LeafCount() - 1, 0);

  tree.Topology()->BinaryIdPreorder(
      [&root_id, &branch_gradient, &height_gradient, leaf_count = tree.LeafCount(),
       &rates = tree.GetRates()](size_t node_id, size_t child0_id, size_t child1_id) {
        if (node_id != root_id) {
          height_gradient[node_id - leaf_count] =
              -branch_gradient[node_id] * rates[node_id];
        }
        if (node_id >= leaf_count) {
          height_gradient[node_id - leaf_count] +=
              branch_gradient[child0_id] * rates[child0_id];
          height_gradient[node_id - leaf_count] +=
              branch_gradient[child1_id] * rates[child1_id];
        }
      });
  return height_gradient;
}

double GetNodePartial(size_t node_id, size_t leaf_count,
                      const std::vector<double> &heights,
                      const std::vector<double> &ratios,
                      const std::vector<double> &bounds) {
  return (heights[node_id] - bounds[node_id]) / ratios[node_id - leaf_count];
}

// Calculate \partial{t_j}/\partial{r_k}
double GetEpochGradientAddition(
    size_t node_id, size_t child_id, size_t leaf_count,
    const std::vector<double> &heights, const std::vector<double> &ratios,
    const std::vector<double> &bounds,
    const std::vector<double> &ratiosGradientUnweightedLogDensity) {
  if (child_id < leaf_count) {
    return 0.0;
  } else if (bounds[node_id] == bounds[child_id]) {
    // child_id and node_id are in the same epoch
    return ratiosGradientUnweightedLogDensity[child_id - leaf_count] *
           ratios[child_id - leaf_count] / ratios[node_id - leaf_count];
  } else {
    // NOT the same epoch
    return ratiosGradientUnweightedLogDensity[child_id - leaf_count] *
           ratios[child_id - leaf_count] / (heights[node_id] - bounds[child_id]) *
           GetNodePartial(node_id, leaf_count, heights, ratios, bounds);
  }
}

std::vector<double> GetLogTimeArray(const RootedTree &tree) {
  size_t leaf_count = tree.LeafCount();
  std::vector<double> log_time(leaf_count - 1, 0);
  const auto &node_bounds = tree.GetNodeBounds();
  const auto &node_heights = tree.GetNodeHeights();
  for (size_t i = 0; i < leaf_count - 2; i++) {
    log_time[i] = 1.0 / (node_heights[leaf_count + i] - node_bounds[leaf_count + i]);
  }
  return log_time;
}

// Update ratio gradient with \partial{t_j}/\partial{r_k}
std::vector<double> UpdateGradientUnWeightedLogDensity(
    const RootedTree &tree, const std::vector<double> &gradient_height) {
  size_t leaf_count = tree.LeafCount();
  size_t root_id = tree.Topology()->Id();
  std::vector<double> ratiosGradientUnweightedLogDensity(leaf_count - 1);
  tree.Topology()->BinaryIdPostorder(
      [&gradient_height, &heights = tree.node_heights_, &ratios = tree.height_ratios_,
       &bounds = tree.GetNodeBounds(), &ratiosGradientUnweightedLogDensity, &leaf_count,
       &root_id](size_t node_id, size_t child0_id, size_t child1_id) {
        if (node_id >= leaf_count && node_id != root_id) {
          ratiosGradientUnweightedLogDensity[node_id - leaf_count] +=
              GetNodePartial(node_id, leaf_count, heights, ratios, bounds) *
              gradient_height[node_id - leaf_count];
          ratiosGradientUnweightedLogDensity[node_id - leaf_count] +=
              GetEpochGradientAddition(node_id, child0_id, leaf_count, heights, ratios,
                                       bounds, ratiosGradientUnweightedLogDensity);
          ratiosGradientUnweightedLogDensity[node_id - leaf_count] +=
              GetEpochGradientAddition(node_id, child1_id, leaf_count, heights, ratios,
                                       bounds, ratiosGradientUnweightedLogDensity);
        }
      });
  return ratiosGradientUnweightedLogDensity;
}

double UpdateHeightParameterGradientUnweightedLogDensity(
    const RootedTree &tree, const std::vector<double> &gradient) {
  size_t leaf_count = tree.LeafCount();
  size_t root_id = tree.Topology()->Id();

  std::vector<double> multiplierArray(leaf_count - 1);
  multiplierArray[root_id - leaf_count] = 1.0;

  tree.Topology()->BinaryIdPreorder(
      [&leaf_count, &ratios = tree.height_ratios_, &multiplierArray](
          size_t node_id, size_t child0_id, size_t child1_id) {
        if (child0_id >= leaf_count) {
          double ratio = ratios[child0_id - leaf_count];
          multiplierArray[child0_id - leaf_count] =
              ratio * multiplierArray[node_id - leaf_count];
        }
        if (child1_id >= leaf_count) {
          double ratio = ratios[child1_id - leaf_count];
          multiplierArray[child1_id - leaf_count] =
              ratio * multiplierArray[node_id - leaf_count];
        }
      });
  double sum = 0.0;
  for (size_t i = 0; i < gradient.size(); i++) {
    sum += gradient[i] * multiplierArray[i];
  }

  return sum;
}

std::vector<double> RatioGradientOfHeightGradient(
    const RootedTree &tree, const std::vector<double> &height_gradient) {
  size_t leaf_count = tree.LeafCount();
  size_t root_id = tree.Topology()->Id();

  // Calculate node ratio gradient
  std::vector<double> gradientLogDensity =
      UpdateGradientUnWeightedLogDensity(tree, height_gradient);

  // Calculate root height gradient
  gradientLogDensity[root_id - leaf_count] =
      UpdateHeightParameterGradientUnweightedLogDensity(tree, height_gradient);

  // Add gradient of log Jacobian determinant
  std::vector<double> log_time = GetLogTimeArray(tree);

  std::vector<double> gradientLogJacobianDeterminant =
      UpdateGradientUnWeightedLogDensity(tree, log_time);
  gradientLogJacobianDeterminant[root_id - leaf_count] =
      UpdateHeightParameterGradientUnweightedLogDensity(tree, log_time);

  for (size_t i = 0; i < gradientLogDensity.size() - 1; i++) {
    gradientLogDensity[i] +=
        gradientLogJacobianDeterminant[i] - 1.0 / tree.height_ratios_[i];
  }

  gradientLogDensity[root_id - leaf_count] +=
      gradientLogJacobianDeterminant[root_id - leaf_count];

  return gradientLogDensity;
}

std::vector<double> RatioGradientOfBranchGradient(
    const RootedTree &tree, const std::vector<double> &branch_gradient) {
  // Calculate node height gradient
  std::vector<double> height_gradient = HeightGradient(tree, branch_gradient);

  return RatioGradientOfHeightGradient(tree, height_gradient);
}

EigenVectorXd RatioGradientOfHeightGradientEigen(
    const RootedTree &tree, EigenConstVectorXdRef height_gradient) {
  std::vector<double> height_gradient_vector(height_gradient.size());
  for (Eigen::Index i = 0; i < height_gradient.size(); ++i) {
    height_gradient_vector[i] = height_gradient(i);
  }
  std::vector<double> vector_output =
      UpdateGradientUnWeightedLogDensity(tree, height_gradient_vector);
  vector_output[vector_output.size() - 1] =
      UpdateHeightParameterGradientUnweightedLogDensity(tree, height_gradient_vector);
  EigenVectorXd eigen_output(vector_output.size());
  for (size_t i = 0; i < vector_output.size(); ++i) {
    eigen_output(i) = vector_output[i];
  }
  return eigen_output;
}
