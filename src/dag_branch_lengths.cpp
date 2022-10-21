// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This class handles DAG General Data and DAG Branch Lengths. Handles storage and
// resizing of data vectors and branch length optimization.

#include "dag_branch_lengths.hpp"

// ** Optimization

void DAGBranchLengths::Optimization(const EdgeId edge_id) {
  Assert(HasDAG(), "DAG must be set to call optimization without only edge_id.");
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  Optimization(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::BrentOptimization(const EdgeId edge_id) {
  Assert(HasDAG(), "DAG must be set to call optimization without only edge_id.");
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  BrentOptimization(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::BrentOptimizationWithGradients(const EdgeId edge_id) {
  Assert(HasDAG(), "DAG must be set to call optimization without only edge_id.");
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  BrentOptimizationWithGradients(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::GradientAscentOptimization(const EdgeId edge_id) {
  Assert(HasDAG(), "DAG must be set to call optimization without only edge_id.");
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  GradientAscentOptimization(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::LogSpaceGradientAscentOptimization(const EdgeId edge_id) {
  Assert(HasDAG(), "DAG must be set to call optimization without only edge_id.");
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  LogSpaceGradientAscentOptimization(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::NewtonOptimization(const EdgeId edge_id) {
  Assert(HasDAG(), "DAG must be set to call optimization without only edge_id.");
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  NewtonOptimization(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::Optimization(const EdgeId edge_id, const NodeId parent_id,
                                    const NodeId child_id) {
  switch (optimization_method_) {
    case Optimization::OptimizationMethod::BrentOptimization:
      return BrentOptimization(edge_id, parent_id, child_id);
    case Optimization::OptimizationMethod::BrentOptimizationWithGradients:
      return BrentOptimizationWithGradients(edge_id, parent_id, child_id);
    case Optimization::OptimizationMethod::GradientAscentOptimization:
      return GradientAscentOptimization(edge_id, parent_id, child_id);
    case Optimization::OptimizationMethod::LogSpaceGradientAscentOptimization:
      return LogSpaceGradientAscentOptimization(edge_id, parent_id, child_id);
    case Optimization::OptimizationMethod::NewtonOptimization:
      return NewtonOptimization(edge_id, parent_id, child_id);
    default:
      Failwith("DAGBranchLengths::Optimization(): Invalid OptimizationMethod given.");
  }
}

void DAGBranchLengths::BrentOptimization(const EdgeId edge_id, const NodeId parent_id,
                                         const NodeId child_id) {
  Assert(neg_llh_f_ != nullptr,
         "NegativeLogLikelihoodFunction must be assigned before calling Brent.");

  // Branch convergence test.
  if (branch_length_differences_(edge_id) < branch_length_difference_threshold_) {
    return;
  }

  auto negative_log_likelihood = [this, edge_id, parent_id,
                                  child_id](double log_branch_length) {
    return neg_llh_f_(edge_id, parent_id, child_id, this->branch_lengths_(edge_id));
  };

  double current_log_branch_length = log(branch_lengths_(edge_id));
  double current_neg_log_likelihood =
      negative_log_likelihood(current_log_branch_length);

  const auto [log_branch_length, neg_log_likelihood] =
      Optimization::BrentMinimize<false>(
          negative_log_likelihood, current_log_branch_length, min_log_branch_length_,
          max_log_branch_length_, significant_digits_for_optimization_,
          max_iter_for_optimization_, step_size_for_log_space_optimization_);

  // Numerical optimization sometimes yields new nllk > current nllk.
  // In this case, we reset the branch length to the previous value.
  if (neg_log_likelihood > current_neg_log_likelihood) {
    branch_lengths_(edge_id) = exp(current_log_branch_length);
  } else {
    branch_lengths_(edge_id) = exp(log_branch_length);
  }
  branch_length_differences_(edge_id) =
      abs(exp(current_log_branch_length) - branch_lengths_(edge_id));
}

void DAGBranchLengths::BrentOptimizationWithGradients(const EdgeId edge_id,
                                                      const NodeId parent_id,
                                                      const NodeId child_id) {
  Assert(neg_llh_f_and_df_ != nullptr,
         "NegativeLogLikelihoodAndDerivativeFunction must be assigned before calling "
         "BrentWithGradients.");

  // Branch convergence test.
  if (branch_length_differences_(edge_id) < branch_length_difference_threshold_) {
    return;
  }

  auto negative_log_likelihood_and_derivative = [this, edge_id, parent_id,
                                                 child_id](double log_branch_length) {
    return this->neg_llh_f_and_df_(edge_id, parent_id, child_id,
                                   this->branch_lengths_(edge_id));
  };

  double current_log_branch_length = log(branch_lengths_(edge_id));
  double current_neg_log_likelihood =
      negative_log_likelihood_and_derivative(current_log_branch_length).first;
  const auto [log_branch_length, neg_log_likelihood] =
      Optimization::BrentMinimize<true>(
          negative_log_likelihood_and_derivative, current_log_branch_length,
          min_log_branch_length_, max_log_branch_length_,
          significant_digits_for_optimization_, max_iter_for_optimization_,
          step_size_for_log_space_optimization_);

  if (neg_log_likelihood > current_neg_log_likelihood) {
    branch_lengths_(edge_id) = exp(current_log_branch_length);
  } else {
    branch_lengths_(edge_id) = exp(log_branch_length);
  }
  branch_length_differences_(edge_id) =
      abs(exp(current_log_branch_length) - branch_lengths_(edge_id));
}

void DAGBranchLengths::GradientAscentOptimization(const EdgeId edge_id,
                                                  const NodeId parent_id,
                                                  const NodeId child_id) {
  Assert(llh_f_and_df_ != nullptr,
         "LogLikelihoodAndDerivativeFunction must be assigned before calling "
         "GradientAscent.");

  auto log_likelihood_and_derivative = [this, edge_id, parent_id,
                                        child_id](double log_branch_length) {
    return this->llh_f_and_df_(edge_id, parent_id, child_id,
                               this->branch_lengths_(edge_id));
  };

  const auto branch_length = Optimization::GradientAscent(
      log_likelihood_and_derivative, branch_lengths_(edge_id),
      significant_digits_for_optimization_, step_size_for_optimization_,
      min_log_branch_length_, max_iter_for_optimization_);
  branch_lengths_(edge_id) = branch_length;
}

void DAGBranchLengths::LogSpaceGradientAscentOptimization(const EdgeId edge_id,
                                                          const NodeId parent_id,
                                                          const NodeId child_id) {
  Assert(llh_f_and_df_ != nullptr,
         "LogLikelihoodAndDerivativeFunction must be assigned before calling "
         "LogSpaceGradientAscent.");

  auto log_likelihood_and_derivative = [this, edge_id, parent_id,
                                        child_id](double log_branch_length) {
    return this->llh_f_and_df_(edge_id, parent_id, child_id,
                               this->branch_lengths_(edge_id));
  };

  const auto branch_length = Optimization::LogSpaceGradientAscent(
      log_likelihood_and_derivative, branch_lengths_(edge_id),
      significant_digits_for_optimization_, step_size_for_log_space_optimization_,
      exp(min_log_branch_length_), max_iter_for_optimization_);
  branch_lengths_(edge_id) = branch_length;
}

void DAGBranchLengths::NewtonOptimization(const EdgeId edge_id, const NodeId parent_id,
                                          const NodeId child_id) {
  Assert(llh_f_and_df_and_ddf_ != nullptr,
         "LogLikelihoodAndDerivativeFunction must be assigned before calling "
         "GradientAscent.");

  if (branch_length_differences_(edge_id) < branch_length_difference_threshold_) {
    return;
  }

  auto log_likelihood_and_first_two_derivatives = [this, edge_id, parent_id,
                                                   child_id](double log_branch_length) {
    return this->llh_f_and_df_and_ddf_(edge_id, parent_id, child_id,
                                       this->branch_lengths_(edge_id));
  };

  double current_log_branch_length = log(branch_lengths_(edge_id));
  const auto log_branch_length = Optimization::NewtonRaphsonOptimization(
      log_likelihood_and_first_two_derivatives, current_log_branch_length,
      significant_digits_for_optimization_, denominator_tolerance_for_newton_,
      min_log_branch_length_, max_log_branch_length_, max_iter_for_optimization_);

  branch_length_differences_(edge_id) =
      abs(exp(current_log_branch_length) - branch_lengths_(edge_id));
}

// ** Optimization Helpers

// void DAGBranchLengths::SetTransitionAndDerivativeMatricesToHaveBranchLength(
//     double branch_length) {
//   diagonal_vector_ = (branch_length * eigenvalues_).array().exp();
//   diagonal_matrix_.diagonal() = diagonal_vector_;
//   transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
//   // Now calculating derivative matrix
//   diagonal_matrix_.diagonal() = eigenvalues_.array() * diagonal_vector_.array();
//   derivative_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
//   // Now calculating hessian matrix
//   diagonal_matrix_.diagonal() =
//       eigenvalues_.array() * eigenvalues_.array() * diagonal_vector_.array();
//   hessian_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
// }

//   DoublePair DAGBranchLengths::LogLikelihoodAndDerivative(
//     const EdgeId edge_id, const NodeId parent_id, const NodeId child_id, const double log_branch_length) {
//   SetTransitionAndDerivativeMatricesToHaveBranchLength(
//       branch_lengths_(edge_id));
//   PreparePerPatternLogLikelihoodsForGPCSP(parent_id.value_, child_id.value_);
//   // The prior is expressed using the current value of q_.
//   // The phylogenetic component of the likelihood is weighted with the number of times
//   // we see the site patterns.
//   const double log_likelihood = per_pattern_log_likelihoods_.dot(site_pattern_weights_);

//   // The per-site likelihood derivative is calculated in the same way as the per-site
//   // likelihood, but using the derivative matrix instead of the transition matrix.
//   // We first prepare two useful vectors _without_ likelihood rescaling, because the
//   // rescalings cancel out in the ratio below.
//   PrepareUnrescaledPerPatternLikelihoodDerivatives(parent_id.value_, child_id.value_);
//   PrepareUnrescaledPerPatternLikelihoods(parent_id.value_, child_id.value_);
//   // If l_i is the per-site likelihood, the derivative of log(l_i) is the derivative
//   // of l_i divided by l_i.
//   per_pattern_likelihood_derivative_ratios_ =
//       per_pattern_likelihood_derivatives_.array() / per_pattern_likelihoods_.array();
//   const double log_likelihood_derivative =
//       per_pattern_likelihood_derivative_ratios_.dot(site_pattern_weights_);
//   return {log_likelihood, log_likelihood_derivative};
// }

// std::tuple<double, double, double> DAGBranchLengths::LogLikelihoodAndFirstTwoDerivatives(
//     const EdgeId edge_id, const NodeId parent_id, const NodeId child_id, const double log_branch_length) {
//   SetTransitionAndDerivativeMatricesToHaveBranchLength(
//       branch_lengths_(edge_id));
//   PreparePerPatternLogLikelihoodsForGPCSP(parent_id.value_, child_id.value_);

//   const double log_likelihood = per_pattern_log_likelihoods_.dot(site_pattern_weights_);

//   // The per-site likelihood derivative is calculated in the same way as the per-site
//   // likelihood, but using the derivative matrix instead of the transition matrix.
//   // We first prepare two useful vectors _without_ likelihood rescaling, because the
//   // rescalings cancel out in the ratio below.
//   PrepareUnrescaledPerPatternLikelihoodDerivatives(parent_id.value_, child_id.value_);
//   PrepareUnrescaledPerPatternLikelihoods(parent_id.value_, child_id.value_);
//   // If l_i is the per-site likelihood, the derivative of log(l_i) is the derivative
//   // of l_i divided by l_i.
//   per_pattern_likelihood_derivative_ratios_ =
//       per_pattern_likelihood_derivatives_.array() / per_pattern_likelihoods_.array();
//   const double log_likelihood_gradient =
//       per_pattern_likelihood_derivative_ratios_.dot(site_pattern_weights_);

//   // Second derivative is calculated the same way, but has an extra term due to
//   // the product rule.

//   PrepareUnrescaledPerPatternLikelihoodSecondDerivatives(parent_id.value_, child_id.value_);

//   per_pattern_likelihood_second_derivative_ratios_ =
//       (per_pattern_likelihood_second_derivatives_.array() *
//            per_pattern_likelihoods_.array() -
//        per_pattern_likelihood_derivatives_.array() *
//            per_pattern_likelihood_derivatives_.array()) /
//       (per_pattern_likelihoods_.array() * per_pattern_likelihoods_.array());

//   const double log_likelihood_hessian =
//       per_pattern_likelihood_second_derivative_ratios_.dot(site_pattern_weights_);

//   return std::make_tuple(log_likelihood, log_likelihood_gradient,
//                          log_likelihood_hessian);
// }
