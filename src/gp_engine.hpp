// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A visitor for GPOperations. See
// https://arne-mertz.de/2018/05/modern-c-features-stdvariant-and-stdvisit/

#ifndef SRC_GP_ENGINE_HPP_
#define SRC_GP_ENGINE_HPP_

#include "eigen_sugar.hpp"
#include "dag_node.hpp"
#include "gp_operation.hpp"
#include "mmapped_plv.hpp"
#include "numerical_utils.hpp"
#include "site_pattern.hpp"
#include "substitution_model.hpp"

class GPEngine {
 public:
  GPEngine(SitePattern site_pattern, size_t pcss_count, std::string mmap_file_path);
  GPEngine(SitePattern site_pattern,
           size_t num_plvs,
           size_t gpcsp_count,
           std::string mmap_file_path);
  
  // These operators mean that we can invoke this class on each of the operations.
  void operator()(const GPOperations::Zero& op);
  void operator()(const GPOperations::SetToStationaryDistribution& op);
  void operator()(const GPOperations::WeightedSumAccumulate& op);
  void operator()(const GPOperations::MarginalLikelihood& op);
  void operator()(const GPOperations::Multiply& op);
  void operator()(const GPOperations::Likelihood& op);
  void operator()(const GPOperations::OptimizeBranchLength& op);
  void operator()(const GPOperations::UpdateSBNProbabilities& op);

  void ProcessOperations(GPOperationVector operations);

  void SetTransitionMatrixToHaveBranchLength(double branch_length);
  void SetTransitionAndDerivativeMatricesToHaveBranchLength(double branch_length);
  void SetTransitionMatrixToHaveBranchLengthAndTranspose(double branch_length);
  const Eigen::Matrix4d& GetTransitionMatrix() { return transition_matrix_; };
  NucleotidePLV GetPLV(size_t idx) { return plvs_[idx]; }
  void PrintPLV(size_t plv_idx);

  void SetBranchLengths(EigenVectorXd branch_lengths) {
    branch_lengths_ = branch_lengths;
  };
  void SetSBNParameters(EigenVectorXd q) {
    q_ = q;
  };
  void ResetLogMarginalLikelihood() {
    log_marginal_likelihood = DOUBLE_NEG_INF;
  }
  double GetLogMarginalLikelihood() {
    return log_marginal_likelihood;
  }
  EigenVectorXd GetBranchLengths() const { return branch_lengths_; };
  EigenVectorXd GetLogLikelihoods() const { return log_likelihoods_; };
  EigenVectorXd GetSBNParameters() const { return q_; };

  DoublePair LogLikelihoodAndDerivative(const GPOperations::OptimizeBranchLength& op);

 private:
  double min_branch_length_ = 1e-6;
  double max_branch_length_ = 3.;
  int significant_digits_for_optimization_ = 6;
  double relative_tolerance_for_optimization_ = 1e-2;
  double step_size_for_optimization_ = 5e-4;
  size_t max_iter_for_optimization_ = 1000;

  double log_marginal_likelihood = DOUBLE_NEG_INF;

  SitePattern site_pattern_;
  size_t plv_count_;
  MmappedNucleotidePLV mmapped_master_plv_;
  // plvs_ store the following:
  // [0, num_nodes): p(s).
  // [num_nodes, 2*num_nodes): phat(s).
  // [2*num_nodes, 3*num_nodes): phat(s_tilde).
  // [3*num_nodes, 4*num_nodes): rhat(s) = rhat(s_tilde).
  // [4*num_nodes, 5*num_nodes): r(s).
  // [5*num_nodes, 6*num_nodes): r(s_tilde).
  NucleotidePLVRefVector plvs_;
  // These parameters are indexed in the same way as sbn_parameters_ in
  // gp_instance.
  EigenVectorXd branch_lengths_;
  EigenVectorXd log_likelihoods_;
  EigenVectorXd q_;

  // Internal "temporaries" useful for likelihood and derivative calculation.
  EigenVectorXd per_pattern_log_likelihoods_;
  EigenVectorXd per_pattern_likelihoods_;
  EigenVectorXd per_pattern_likelihood_derivatives_;
  EigenVectorXd per_pattern_likelihood_derivative_ratios_;

  // When we change from JC69Model, check that we are actually doing transpose in
  // leafward calculations.
  JC69Model substitution_model_;
  Eigen::Matrix4d eigenmatrix_ = substitution_model_.GetEigenvectors().reshaped(4, 4);
  Eigen::Matrix4d inverse_eigenmatrix_ =
      substitution_model_.GetInverseEigenvectors().reshaped(4, 4);
  Eigen::Vector4d eigenvalues_ = substitution_model_.GetEigenvalues();
  Eigen::DiagonalMatrix<double, 4> diagonal_matrix_;
  Eigen::Matrix4d transition_matrix_;
  Eigen::Matrix4d derivative_matrix_;
  Eigen::Vector4d stationary_distribution_ = substitution_model_.GetFrequencies();
  EigenVectorXd site_pattern_weights_;

  void InitializePLVsWithSitePatterns();
  void BrentOptimization(const GPOperations::OptimizeBranchLength& op);
  void GradientAscentOptimization(const GPOperations::OptimizeBranchLength& op);

//  inline double LogLikelihood(size_t src1_idx, size_t src2_idx) {
//    per_pattern_log_likelihoods_ =
//        (plvs_.at(src1_idx).transpose() * plvs_.at(src2_idx)).diagonal().array().log();
//    return per_pattern_log_likelihoods_.dot(site_pattern_weights_);
//  }

  inline void PreparePerPatternLikelihoodDerivatives(size_t src1_idx, size_t src2_idx) {
    per_pattern_likelihood_derivatives_ =
        (plvs_.at(src1_idx).transpose() * derivative_matrix_ * plvs_.at(src2_idx)).diagonal().array();
  }

  inline void PreparePerPatternLikelihoods(size_t src1_idx, size_t src2_idx) {
    per_pattern_likelihoods_ =
        (plvs_.at(src1_idx).transpose() * transition_matrix_ * plvs_.at(src2_idx)).diagonal().array();
  }

  inline DoublePair LogLikelihoodAndDerivativeFromPreparations() {
    per_pattern_log_likelihoods_ = per_pattern_likelihoods_.array().log();
    double log_likelihood = per_pattern_log_likelihoods_.dot(site_pattern_weights_);
    // If l_i is the per-site likelihood, the derivative of log(l_i) is the derivative
    // of l_i divided by l_i.
    per_pattern_likelihood_derivative_ratios_ =
        per_pattern_likelihood_derivatives_.array() / per_pattern_likelihoods_.array();
    // We weight this with the number of times we see the site patterns as for the log
    // likelihood.
    double log_likelihood_derivative =
        per_pattern_likelihood_derivative_ratios_.dot(site_pattern_weights_);
    return {log_likelihood, log_likelihood_derivative};
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("GPEngine") {
  SitePattern hello_site_pattern = SitePattern::HelloSitePattern();
  GPEngine engine(hello_site_pattern, 5, "_ignore/mmapped_plv.data");
  engine.SetTransitionMatrixToHaveBranchLength(0.75);
  // Computed directly:
  // https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_%28Jukes_and_Cantor_1969%29
  CHECK(fabs(0.52590958087 - engine.GetTransitionMatrix()(0, 0)) < 1e-10);
  CHECK(fabs(0.1580301397 - engine.GetTransitionMatrix()(0, 1)) < 1e-10);
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_GP_ENGINE_HPP_
