// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Calculation of the ratio and root height gradient, adapted from BEAST.
// https://github.com/beast-dev/beast-mcmc
// Credit to Xiang Ji and Marc Suchard.
//
// Because this is code adapted from elsewhere, at least for the time being the naming
// conventions are a little different: pascalCase is allowed for variables.

#include "rooted_gradient_transforms.hpp"

#include <numeric>

#include "rooted_tree.hpp"

// \partial{L}/\partial{t_k} = \sum_j \partial{L}/\partial{b_j}
// \partial{b_j}/\partial{t_k}
std::vector<double> RootedGradientTransforms::HeightGradient(
    const RootedTree &tree, const std::vector<double> &branch_gradient) {
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

double RootedGradientTransforms::GetNodePartial(size_t node_id, size_t leaf_count,
                                                const std::vector<double> &heights,
                                                const std::vector<double> &ratios,
                                                const std::vector<double> &bounds) {
  return (heights[node_id] - bounds[node_id]) / ratios[node_id - leaf_count];
}

// Calculate \partial{t_j}/\partial{r_k}
double RootedGradientTransforms::GetEpochGradientAddition(
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

std::vector<double> RootedGradientTransforms::GetLogTimeArray(const RootedTree &tree) {
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
std::vector<double> RootedGradientTransforms::UpdateGradientUnWeightedLogDensity(
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

double RootedGradientTransforms::UpdateHeightParameterGradientUnweightedLogDensity(
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

std::vector<double> RootedGradientTransforms::GradientLogDeterminantJacobian(
    const RootedTree &tree) {
  size_t leaf_count = tree.LeafCount();
  size_t root_id = tree.Topology()->Id();

  std::vector<double> log_time = GetLogTimeArray(tree);

  std::vector<double> gradient_log_jacobian_determinant =
      UpdateGradientUnWeightedLogDensity(tree, log_time);

  gradient_log_jacobian_determinant[root_id - leaf_count] =
      UpdateHeightParameterGradientUnweightedLogDensity(tree, log_time);

  for (size_t i = 0; i < gradient_log_jacobian_determinant.size() - 1; i++) {
    gradient_log_jacobian_determinant[i] -= 1.0 / tree.height_ratios_[i];
  }

  return gradient_log_jacobian_determinant;
}

std::vector<double> RootedGradientTransforms::RatioGradientOfHeightGradient(
    const RootedTree &tree, const std::vector<double> &height_gradient) {
  size_t leaf_count = tree.LeafCount();
  size_t root_id = tree.Topology()->Id();

  // Calculate node ratio gradient
  std::vector<double> gradient_log_density =
      UpdateGradientUnWeightedLogDensity(tree, height_gradient);

  // Calculate root height gradient
  gradient_log_density[root_id - leaf_count] =
      UpdateHeightParameterGradientUnweightedLogDensity(tree, height_gradient);

  return gradient_log_density;
}

std::vector<double> RootedGradientTransforms::RatioGradientOfBranchGradient(
    const RootedTree &tree, const std::vector<double> &branch_gradient) {
  size_t leaf_count = tree.LeafCount();
  size_t root_id = tree.Topology()->Id();

  // Calculate node height gradient
  std::vector<double> height_gradient = HeightGradient(tree, branch_gradient);

  // Calculate ratios and root height gradient
  std::vector<double> gradient_log_density =
      RatioGradientOfHeightGradient(tree, height_gradient);

  // Calculate gradient of log Jacobian determinant
  std::vector<double> gradient_log_jacobian_determinant =
      GradientLogDeterminantJacobian(tree);

  for (size_t i = 0; i < gradient_log_jacobian_determinant.size() - 1; i++) {
    gradient_log_density[i] += gradient_log_jacobian_determinant[i];
  }

  gradient_log_density[root_id - leaf_count] +=
      gradient_log_jacobian_determinant[root_id - leaf_count];

  return gradient_log_density;
}

std::vector<double> RootedGradientTransforms::RatioGradientOfBranchGradient(
    const RootedTree &tree, const std::vector<double> &branch_gradient,
    const std::optional<PhyloFlags> flags) {
  size_t leaf_count = tree.LeafCount();
  size_t root_id = tree.Topology()->Id();

  // Calculate node height gradient
  std::vector<double> height_gradient = HeightGradient(tree, branch_gradient);

  // Calculate ratios and root height gradient
  std::vector<double> gradient_log_density =
      RatioGradientOfHeightGradient(tree, height_gradient);

  if (PhyloFlags::IsFlagSet(
          flags, PhyloGradientFlagOptions::include_log_det_jacobian_gradient_)) {
    std::vector<double> gradient_log_jacobian_determinant =
        GradientLogDeterminantJacobian(tree);

    for (size_t i = 0; i < gradient_log_jacobian_determinant.size() - 1; i++) {
      gradient_log_density[i] += gradient_log_jacobian_determinant[i];
    }

    gradient_log_density[root_id - leaf_count] +=
        gradient_log_jacobian_determinant[root_id - leaf_count];
  }

  return gradient_log_density;
}

EigenVectorXd RootedGradientTransforms::RatioGradientOfHeightGradientEigen(
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

double RootedGradientTransforms::LogDetJacobianHeightTransform(const RootedTree &tree) {
  double log_det_jacobian = 0.0;
  size_t leaf_count = tree.LeafCount();
  tree.Topology()->TripleIdPreorderBifurcating(
      [&log_det_jacobian, &tree, leaf_count](int node_id, int sister_id,
                                             int parent_id) {
        if (size_t(node_id) >=
            leaf_count) {  // Only add to computation if node is not a leaf.
          // Account for the jacobian of this branch's height transform.
          log_det_jacobian +=
              std::log(tree.node_heights_[parent_id] - tree.node_bounds_[node_id]);
        }
      });
  return log_det_jacobian;
}
