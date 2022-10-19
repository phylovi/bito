// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This class handles DAG General Data and DAG Branch Lengths. Handles storage and
// resizing of data vectors and branch length optimization.

#include "dag_branch_lengths.hpp"

void DAGBranchLengths::Optimization(const EdgeId edge_id) {
  switch (optimization_method_) {
    case Optimization::OptimizationMethod::BrentOptimization:
      return BrentOptimization(edge_id);
    case Optimization::OptimizationMethod::BrentOptimizationWithGradients:
      return BrentOptimizationWithGradients(edge_id);
    case Optimization::OptimizationMethod::GradientAscentOptimization:
      return GradientAscentOptimization(edge_id);
    case Optimization::OptimizationMethod::LogSpaceGradientAscentOptimization:
      return LogSpaceGradientAscentOptimization(edge_id);
    case Optimization::OptimizationMethod::NewtonOptimization:
      return NewtonOptimization(edge_id);
    default:
      Failwith("DAGBranchLengths::Optimization(): Invalid OptimizationMethod given.");
  }
}

void DAGBranchLengths::BrentOptimization(const EdgeId edge_id) {
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  BrentOptimization(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::BrentOptimizationWithGradients(const EdgeId edge_id) {
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  BrentOptimizationWithGradients(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::GradientAscentOptimization(const EdgeId edge_id) {
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  GradientAscentOptimization(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::LogSpaceGradientAscentOptimization(const EdgeId edge_id) {
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  LogSpaceGradientAscentOptimization(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::NewtonOptimization(const EdgeId edge_id) {
  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  NewtonOptimization(edge.GetId(), edge.GetParent(), edge.GetChild());
}

void DAGBranchLengths::BrentOptimization(const EdgeId edge_id, const NodeId parent_id,
                                         const NodeId child_id) {
  Assert(neg_llh_f_ != nullptr,
         "NegativeLogLikelihoodFunction must be assigned before calling Brent.");

  // Branch convergence test.
  if (branch_length_differences_(edge_id) < branch_length_difference_threshold_) {
    return;
  }

  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  auto negative_log_likelihood = [this, edge](double log_branch_length) {
    return neg_llh_f_(edge.GetId(), edge.GetParent(), edge.GetParent(),
                      this->branch_lengths_(edge.GetId()));
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

  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  auto negative_log_likelihood_and_derivative = [this, edge](double log_branch_length) {
    return this->neg_llh_f_and_df_(edge.GetId(), edge.GetParent(), edge.GetChild(),
                                   this->branch_lengths_(edge.GetId()));
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

  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  auto log_likelihood_and_derivative = [this, edge](double log_branch_length) {
    return this->llh_f_and_df_(edge.GetId(), edge.GetParent(), edge.GetChild(),
                               this->branch_lengths_(edge.GetId()));
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

  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  auto log_likelihood_and_derivative = [this, edge](double log_branch_length) {
    return this->llh_f_and_df_(edge.GetId(), edge.GetParent(), edge.GetChild(),
                               this->branch_lengths_(edge.GetId()));
  };

  const auto branch_length = Optimization::LogSpaceGradientAscent(
      log_likelihood_and_derivative, branch_lengths_(edge_id),
      significant_digits_for_optimization_, step_size_for_log_space_optimization_,
      exp(min_log_branch_length_), max_iter_for_optimization_);
  branch_lengths_(edge_id) = branch_length;
}

void DAGBranchLengths::NewtonOptimization(const EdgeId edge_id, const NodeId parent_id,
                                          const NodeId child_id) {
  Assert(llh_f_and_df_ != nullptr,
         "LogLikelihoodAndDerivativeFunction must be assigned before calling "
         "GradientAscent.");

  if (branch_length_differences_(edge_id) < branch_length_difference_threshold_) {
    return;
  }

  const auto& edge = GetDAG().GetDAGEdge(edge_id);
  auto log_likelihood_and_first_two_derivatives = [this,
                                                   edge](double log_branch_length) {
    return this->llh_f_and_df_and_ddf_(edge.GetId(), edge.GetParent(), edge.GetChild(),
                                       this->branch_lengths_(edge.GetId()));
  };

  double current_log_branch_length = log(branch_lengths_(edge_id));
  const auto log_branch_length = Optimization::NewtonRaphsonOptimization(
      log_likelihood_and_first_two_derivatives, current_log_branch_length,
      significant_digits_for_optimization_, denominator_tolerance_for_newton_,
      min_log_branch_length_, max_log_branch_length_, max_iter_for_optimization_);

  branch_lengths_(edge_id) = exp(log_branch_length);

  branch_length_differences_(edge_id) =
      abs(exp(current_log_branch_length) - branch_lengths_(edge_id));
}
