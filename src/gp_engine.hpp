// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A visitor for GPOperations. See
// https://arne-mertz.de/2018/05/modern-c-features-stdvariant-and-stdvisit/

#ifndef SRC_GP_ENGINE_HPP_
#define SRC_GP_ENGINE_HPP_

#include "eigen_sugar.hpp"
#include "gp_operation.hpp"
#include "mmapped_plv.hpp"
#include "numerical_utils.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "site_pattern.hpp"
#include "substitution_model.hpp"

class GPEngine {
 public:
  GPEngine(SitePattern site_pattern, size_t plv_count, size_t gpcsp_count,
           const std::string& mmap_file_path, double rescaling_threshold);

  // These operators mean that we can invoke this class on each of the operations.
  void operator()(const GPOperations::Zero& op);
  void operator()(const GPOperations::SetToStationaryDistribution& op);
  void operator()(const GPOperations::IncrementWithWeightedEvolvedPLV& op);
  void operator()(const GPOperations::IncrementMarginalLikelihood& op);
  void operator()(const GPOperations::Multiply& op);
  void operator()(const GPOperations::Likelihood& op);
  void operator()(const GPOperations::OptimizeBranchLength& op);
  void operator()(const GPOperations::UpdateSBNProbabilities& op);
  void operator()(const GPOperations::PrepForMarginalization& op);

  void ProcessOperations(GPOperationVector operations);

  void SetTransitionMatrixToHaveBranchLength(double branch_length);
  void SetTransitionAndDerivativeMatricesToHaveBranchLength(double branch_length);
  void SetTransitionMatrixToHaveBranchLengthAndTranspose(double branch_length);
  const Eigen::Matrix4d& GetTransitionMatrix() { return transition_matrix_; };
  void PrintPLV(size_t plv_idx);

  void SetBranchLengths(EigenVectorXd branch_lengths) {
    branch_lengths_ = std::move(branch_lengths);
  };
  void SetBranchLengthsToConstant(double branch_length) {
    branch_lengths_.setConstant(branch_length);
  };
  void SetSBNParameters(EigenVectorXd q) { q_ = std::move(q); };
  void ResetLogMarginalLikelihood() {
    log_marginal_likelihood_.setConstant(DOUBLE_NEG_INF);
  }
  double GetLogMarginalLikelihood() const {
    return (log_marginal_likelihood_.array() * site_pattern_weights_.array()).sum();
  }
  EigenVectorXd GetBranchLengths() const { return branch_lengths_; };
  // EigenVectorXd GetLogLikelihoods() const { return log_likelihoods_; };
  EigenVectorXd GetSBNParameters() const { return q_; };

  // Use branch lengths from loaded sample as a starting point for optimization.
  void HotStartBranchLengths(const RootedTreeCollection& tree_collection,
                             const BitsetSizeMap& indexer);

  DoublePair LogLikelihoodAndDerivative(const GPOperations::OptimizeBranchLength& op);

  static constexpr double default_rescaling_threshold_ = 1e-40;
  static constexpr double default_branch_length_ = 0.1;

  double PLVByteCount() const { return mmapped_master_plv_.ByteCount(); };

 private:
  static constexpr double min_branch_length_ = 1e-6;
  static constexpr double max_branch_length_ = 3.;

  int significant_digits_for_optimization_ = 6;
  double relative_tolerance_for_optimization_ = 1e-2;
  double step_size_for_optimization_ = 5e-4;
  size_t max_iter_for_optimization_ = 1000;

  EigenVectorXd log_marginal_likelihood_;

  SitePattern site_pattern_;
  size_t plv_count_;
  const double rescaling_threshold_;
  const double log_rescaling_threshold_;
  MmappedNucleotidePLV mmapped_master_plv_;
  // plvs_ store the following (see GPDAG::GetPLVIndexStatic):
  // [0, num_nodes): p(s).
  // [num_nodes, 2*num_nodes): phat(s).
  // [2*num_nodes, 3*num_nodes): phat(s_tilde).
  // [3*num_nodes, 4*num_nodes): rhat(s) = rhat(s_tilde).
  // [4*num_nodes, 5*num_nodes): r(s).
  // [5*num_nodes, 6*num_nodes): r(s_tilde).
  NucleotidePLVRefVector plvs_;
  EigenVectorXi rescaling_counts_;
  // These parameters are indexed in the same way as sbn_parameters_ in
  // gp_instance.
  EigenVectorXd branch_lengths_;
  EigenMatrixXd log_likelihoods_;
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
  Eigen::Vector4d diagonal_vector_;
  Eigen::DiagonalMatrix<double, 4> diagonal_matrix_;
  Eigen::Matrix4d transition_matrix_;
  Eigen::Matrix4d derivative_matrix_;
  Eigen::Vector4d stationary_distribution_ = substitution_model_.GetFrequencies();
  EigenVectorXd site_pattern_weights_;

  void InitializePLVsWithSitePatterns();

  void RescalePLV(size_t plv_idx, int amount);
  void AssertPLVIsFinite(size_t plv_idx, const std::string& message) const;
  std::pair<double, double> PLVMinMax(size_t plv_idx) const;
  // If a PLV all entries smaller than rescaling_threshold_ then rescale it up and
  // increment the corresponding entry in rescaling_counts_.
  void RescalePLVIfNeeded(size_t plv_idx);
  double LogRescalingFor(size_t plv_idx);

  void BrentOptimization(const GPOperations::OptimizeBranchLength& op);
  void GradientAscentOptimization(const GPOperations::OptimizeBranchLength& op);

  inline void PrepareUnrescaledPerPatternLikelihoodDerivatives(size_t src1_idx,
                                                               size_t src2_idx) {
    per_pattern_likelihood_derivatives_ =
        (plvs_.at(src1_idx).transpose() * derivative_matrix_ * plvs_.at(src2_idx))
            .diagonal()
            .array();
  }

  inline void PrepareUnrescaledPerPatternLikelihoods(size_t src1_idx, size_t src2_idx) {
    per_pattern_likelihoods_ =
        (plvs_.at(src1_idx).transpose() * transition_matrix_ * plvs_.at(src2_idx))
            .diagonal()
            .array();
  }

  inline void PreparePerPatternLogLikelihoods(size_t src1_idx, size_t src2_idx) {
    per_pattern_log_likelihoods_ =
        (plvs_.at(src1_idx).transpose() * transition_matrix_ * plvs_.at(src2_idx))
            .diagonal()
            .array()
            .log() +
        LogRescalingFor(src1_idx) + LogRescalingFor(src2_idx);
  }

  inline void PreparePerPatternLogLikelihoods(size_t gpcsp_idx, size_t src1_idx,
                                              size_t src2_idx) {
    per_pattern_log_likelihoods_ = (q_[gpcsp_idx] * plvs_.at(src1_idx).transpose() *
                                    transition_matrix_ * plvs_.at(src2_idx))
                                       .diagonal()
                                       .array()
                                       .log() +
                                   LogRescalingFor(src1_idx) +
                                   LogRescalingFor(src2_idx);
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("GPEngine") {
  SitePattern hello_site_pattern = SitePattern::HelloSitePattern();
  GPEngine engine(hello_site_pattern, 6 * 5, 5, "_ignore/mmapped_plv.data",
                  GPEngine::default_rescaling_threshold_);
  engine.SetTransitionMatrixToHaveBranchLength(0.75);
  // Computed directly:
  // https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_%28Jukes_and_Cantor_1969%29
  CHECK(fabs(0.52590958087 - engine.GetTransitionMatrix()(0, 0)) < 1e-10);
  CHECK(fabs(0.1580301397 - engine.GetTransitionMatrix()(0, 1)) < 1e-10);
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_GP_ENGINE_HPP_
