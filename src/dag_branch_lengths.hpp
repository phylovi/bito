// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This class handles DAG General Data and DAG Branch Lengths. Handles storage and
// resizing of data vectors and branch length optimization.

#pragma once

#include "gp_dag.hpp"
#include "dag_data.hpp"
#include "optimization.hpp"
#include "substitution_model.hpp"

class DAGBranchLengths {
 public:
  // LogLikelihood functions take the following args: (1) focal edge id, (2) rootward
  // node id, (3) leafward node id, (4) log branch length of focal edge.
  using LogLikelihoodAndDerivativeFunc =
      std::function<DoublePair(EdgeId, NodeId, NodeId, double)>;
  using LogLikelihoodAndFirstTwoDerivativesFunc =
      std::function<std::tuple<double, double, double>(EdgeId, NodeId, NodeId, double)>;
  using NegLogLikelihoodFunc = std::function<double(EdgeId, NodeId, NodeId, double)>;
  using NegLogLikelihoodAndDerivativeFunc =
      std::function<DoublePair(EdgeId, NodeId, NodeId, double)>;

  DAGBranchLengths(
      const size_t count,
      std::optional<Optimization::OptimizationMethod> method = std::nullopt)
      : branch_lengths_(count), branch_length_differences_(count) {
    if (method.has_value()) {
      SetOptimizationMethod(method.value());
    }
    branch_lengths_.SetDefaultValue(default_branch_length_);
    branch_length_differences_.SetDefaultValue(default_branch_length_);
  }
  DAGBranchLengths(
      GPDAG& dag, std::optional<Optimization::OptimizationMethod> method = std::nullopt)
      : branch_lengths_(dag), branch_length_differences_(dag) {
    if (method.has_value()) {
      SetOptimizationMethod(method.value());
    }
    branch_lengths_.SetDefaultValue(default_branch_length_);
    branch_length_differences_.SetDefaultValue(default_branch_length_);
  }

  // ** Counts

  size_t GetCount() { return GetBranchLengths().GetCount(); }
  size_t GetSpareCount() { return GetBranchLengths().GetSpareCount(); }
  size_t GetAllocCount() { return GetBranchLengths().GetAllocCount(); }

  void SetCount(const size_t count) {
    GetBranchLengths().SetCount(count);
    GetBranchDifferences().SetCount(count);
  }
  void SetSpareCount(const size_t count) {
    GetBranchLengths().SetSpareCount(count);
    GetBranchDifferences().SetSpareCount(count);
  }
  void SetAllocCount(const size_t count) {
    GetBranchLengths().SetAllocCount(count);
    GetBranchDifferences().SetAllocCount(count);
  }

  // ** Access

  // Get the size of the data vector.
  size_t size() { return branch_lengths_.GetData().size(); }

  // Get full data vectors.
  DAGEdgeDoubleData& GetBranchLengths() { return branch_lengths_; }
  const DAGEdgeDoubleData& GetBranchLengths() const { return branch_lengths_; }
  DAGEdgeDoubleData& GetBranchDifferences() { return branch_length_differences_; }
  const DAGEdgeDoubleData& GetBranchDifferences() const {
    return branch_length_differences_;
  }

  // Access elements in data vector.
  double& Get(const EdgeId edge_id) { return branch_lengths_(edge_id); }
  double& operator()(const EdgeId edge_id) { return Get(edge_id); }
  double& GetBranchDifference(const EdgeId edge_id) {
    return branch_length_differences_(edge_id);
  }

  // Check if there is a reference vector.
  const bool HasDAG() const { return branch_lengths_.HasDAG(); }
  const GPDAG& GetDAG() const { return branch_lengths_.GetDAG(); }

  void SetOptimizationMethod(const Optimization::OptimizationMethod method) {
    optimization_method_ = method;
  }
  void SetSignificantDigitsForOptimization(int significant_digits) {
    significant_digits_for_optimization_ = significant_digits;
  }

  void SetLogLikelihoodAndDerivativeFunc(LogLikelihoodAndDerivativeFunc llh_f_and_df) {
    llh_f_and_df_ = llh_f_and_df;
  }
  void SetLogLikelihoodAndFirstTwoDerivativesFunc(
      LogLikelihoodAndFirstTwoDerivativesFunc llh_f_and_df_and_ddf) {
    llh_f_and_df_and_ddf_ = llh_f_and_df_and_ddf;
  }
  void SetNegLogLikelihoodFunc(NegLogLikelihoodFunc neg_llh_f) {
    neg_llh_f_ = neg_llh_f;
  }
  void SetNegLogLikelihoodAndDerivativeFunc(
      NegLogLikelihoodAndDerivativeFunc neg_llh_f_and_df) {
    neg_llh_f_and_df_ = neg_llh_f_and_df;
  }

  Optimization::OptimizationMethod GetOptimizationMethod() {
    return optimization_method_;
  }

  // ** Resize

  void Resize(const size_t count,
              std::optional<const Reindexer> reindexer = std::nullopt,
              std::optional<const size_t> explicit_alloc = std::nullopt) {
    branch_lengths_.Resize(count, reindexer, explicit_alloc);
    branch_length_differences_.Resize(count, reindexer, explicit_alloc);
  }

  void Resize(std::optional<const Reindexer> reindexer = std::nullopt,
              std::optional<const size_t> explicit_alloc = std::nullopt) {
    branch_lengths_.Resize(reindexer, explicit_alloc);
    branch_length_differences_.Resize(reindexer, explicit_alloc);
  }

  // ** Optimization

  // These optimization functions can be called without an assigned reference DAG.
  void Optimization(const EdgeId edge_id, const NodeId parent_id,
                    const NodeId child_id);
  void BrentOptimization(const EdgeId edge_id, const NodeId parent_id,
                         const NodeId child_id);
  void BrentOptimizationWithGradients(const EdgeId edge_id, const NodeId parent_id,
                                      const NodeId child_id);
  void GradientAscentOptimization(const EdgeId edge_id, const NodeId parent_id,
                                  const NodeId child_id);
  void LogSpaceGradientAscentOptimization(const EdgeId edge_id, const NodeId parent_id,
                                          const NodeId child_id);
  void NewtonOptimization(const EdgeId edge_id, const NodeId parent_id,
                          const NodeId child_id);

  // These optimization functions require an assigned reference DAG.
  void Optimization(const EdgeId edge_id);
  void BrentOptimization(const EdgeId edge_id);
  void BrentOptimizationWithGradients(const EdgeId edge_id);
  void GradientAscentOptimization(const EdgeId edge_id);
  void LogSpaceGradientAscentOptimization(const EdgeId edge_id);
  void NewtonOptimization(const EdgeId edge_id);

  // ** Optimization Helpers

void SetTransitionAndDerivativeMatricesToHaveBranchLength(
    double branch_length);
  DoublePair LogLikelihoodAndDerivative(
    EdgeId edge_id, NodeId parent_id, NodeId child_id, double log_branch_length);
std::tuple<double, double, double> LogLikelihoodAndFirstTwoDerivatives(
    EdgeId edge_id, NodeId parent_id, NodeId child_id, double log_branch_length);

 protected:
  // Branch lengths.
  DAGEdgeDoubleData branch_lengths_;
  // Tracks differences between branch length over different iterations of
  // optimization to test for convergence.
  DAGEdgeDoubleData branch_length_differences_;
  // Value to initialize new branch lengths to.
  double default_branch_length_ = 0.1;

  // Method used in branch length optimization.
  Optimization::OptimizationMethod optimization_method_;
  // Absolute lower bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double min_log_branch_length_ = -13.9;
  // Absolute upper bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double max_log_branch_length_ = 1.1;
  // Precision used for checking convergence of branch length optimization.
  // In the non-Brent optimization methods, significant digits will be used to
  // determine convergence through relative tolerance, i.e. measuring difference
  // from previous branch length values until the absolute difference is below
  // 10^(-significant_digits_for_optimization_).
  // Brent optimization does not define convergence through relative tolerance,
  // rather convergence based on tightness of the brackets that it adapts during
  // optimization. This variable thus represents the "number of bits precision to
  // which the minima should be found". When testing on sample datasets, we found that
  // setting the value to 10 was a good compromise between speed and precision for
  // Brent. See more on Brent optimization here:
  // https://www.boost.org/doc/libs/1_79_0/libs/math/doc/html/math_toolkit/brent_minima.html
  int significant_digits_for_optimization_ = 10;
  double relative_tolerance_for_optimization_ = 1e-4;
  double denominator_tolerance_for_newton_ = 1e-10;
  double step_size_for_optimization_ = 5e-4;
  double step_size_for_log_space_optimization_ = 1.0005;
  // Number of iterations allowed for branch length optimization.
  size_t max_iter_for_optimization_ = 1000;
  double branch_length_difference_threshold_ = 1e-15;

  // ** Model

  // When we change from JC69Model, check that we are actually doing transpose in
  // leafward calculations.
  // JC69Model substitution_model_;
  // Eigen::Matrix4d eigenmatrix_ = substitution_model_.GetEigenvectors().reshaped(4, 4);
  // Eigen::Matrix4d inverse_eigenmatrix_ =
  //     substitution_model_.GetInverseEigenvectors().reshaped(4, 4);
  // Eigen::Vector4d eigenvalues_ = substitution_model_.GetEigenvalues();
  // Eigen::Vector4d diagonal_vector_;
  // Eigen::DiagonalMatrix<double, 4> diagonal_matrix_;
  // Eigen::Matrix4d transition_matrix_;
  // Eigen::Matrix4d derivative_matrix_;
  // Eigen::Matrix4d hessian_matrix_;
  // Eigen::Vector4d stationary_distribution_ = substitution_model_.GetFrequencies();
  // EigenVectorXd site_pattern_weights_;

  // Optimization Helper functions
  LogLikelihoodAndDerivativeFunc llh_f_and_df_ = nullptr;
  LogLikelihoodAndFirstTwoDerivativesFunc llh_f_and_df_and_ddf_ = nullptr;
  NegLogLikelihoodFunc neg_llh_f_ = nullptr;
  NegLogLikelihoodAndDerivativeFunc neg_llh_f_and_df_ = nullptr;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("DAG Data") { std::cout << " === DAG Data TEST_CASE ===" << std::endl; }

#endif  // DOCTEST_LIBRARY_INCLUDED
