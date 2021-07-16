// Copyright 2019-2021 libsbn project contributors.
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
#include "quartet_hybrid_request.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "site_pattern.hpp"
#include "substitution_model.hpp"

class GPEngine {
 public:
  GPEngine(SitePattern site_pattern, size_t plv_count, size_t gpcsp_count,
           const std::string& mmap_file_path, double rescaling_threshold,
           EigenVectorXd sbn_prior, EigenVectorXd unconditional_node_probabilities,
           EigenVectorXd inverted_sbn_prior);

  // These operators mean that we can invoke this class on each of the operations.
  void operator()(const GPOperations::ZeroPLV& op);
  void operator()(const GPOperations::SetToStationaryDistribution& op);
  void operator()(const GPOperations::IncrementWithWeightedEvolvedPLV& op);
  void operator()(const GPOperations::ResetMarginalLikelihood& op);
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

  void SetBranchLengths(EigenVectorXd branch_lengths);
  void SetBranchLengthsToConstant(double branch_length);
  void ResetLogMarginalLikelihood();
  double GetLogMarginalLikelihood() const;
  EigenVectorXd GetBranchLengths() const;
  // This function returns a vector indexed by GPCSP such that the i-th entry
  // stores the log of the across-sites product of
  // (the marginal likelihood conditioned on a given GPCSP) *
  //     (the unconditional probability of i's parent subsplit).
  // That is, it's sum_m r^m(t) P(t -> s) p^m(s).
  // See lem:PerPCSPMarginalLikelihood.
  // #288 rename?
  EigenVectorXd GetPerGPCSPLogLikelihoods() const;
  // This override of GetPerGPCSPLogLikelihoods computes the marginal log
  // likelihood for GPCSPs in the range [start, start + length).
  EigenVectorXd GetPerGPCSPLogLikelihoods(size_t start, size_t length) const;
  // This is the full marginal likelihood sum restricted to trees containing a PCSP.
  // When we sum the log of eq:PerGPCSPComponentsOfFullMarginal over the sites, we get
  // out a term that is the number of sites times the log of the prior conditional PCSP
  // probability.
  EigenVectorXd GetPerGPCSPComponentsOfFullLogMarginal() const;
  // #288 reconsider this name
  EigenConstMatrixXdRef GetLogLikelihoodMatrix() const;
  EigenConstVectorXdRef GetHybridMarginals() const;
  EigenConstVectorXdRef GetSBNParameters() const;

  // Calculate a vector of likelihoods, one for each summand of the hybrid marginal.
  EigenVectorXd CalculateQuartetHybridLikelihoods(const QuartetHybridRequest& request);
  // Calculate the actual hybrid marginal and store it in the corresponding entry of
  // hybrid_marginal_log_likelihoods_.
  void ProcessQuartetHybridRequest(const QuartetHybridRequest& request);

  void PrintPLV(size_t plv_idx);

  // Use branch lengths from loaded sample as a starting point for optimization.
  void HotStartBranchLengths(const RootedTreeCollection& tree_collection,
                             const BitsetSizeMap& indexer);

  DoublePair LogLikelihoodAndDerivative(const GPOperations::OptimizeBranchLength& op);
  std::tuple<double, double, double> LogLikelihoodAndFirstTwoDerivatives(
      const GPOperations::OptimizeBranchLength& op);

  static constexpr double default_rescaling_threshold_ = 1e-40;
  static constexpr double default_branch_length_ = 0.1;

  double PLVByteCount() const { return mmapped_master_plv_.ByteCount(); };

 private:
  static constexpr double min_log_branch_length_ = -13.9;
  static constexpr double max_log_branch_length_ = 1.1;

  int significant_digits_for_optimization_ = 6;
  double relative_tolerance_for_optimization_ = 1e-3;
  double step_size_for_optimization_ = 5e-4;
  double step_size_for_log_space_optimization_ = 1.0003;
  size_t max_iter_for_optimization_ = 10000;

  int montonicity_const_for_adaptive_stepsize_ = 10;
  double adaptive_stepsize_uniformity_bound_ = 1e-10;
  double threshold_const_for_adaptive_stepsize_ = 1e-4;
  double rescaling_const_for_adaptive_stepsize_ = 0.5;

  // The length of this vector is equal to the number of site patterns.
  // Entry j stores the marginal log likelihood over all trees at site pattern
  // j.
  EigenVectorXd log_marginal_likelihood_;

  // This vector is indexed by the GPCSPs and stores the hybrid marginals if they are
  // available.
  EigenVectorXd hybrid_marginal_log_likelihoods_;

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
  // branch_lengths_, q_, etc. are indexed in the same way as sbn_parameters_ in
  // gp_instance.
  EigenVectorXd branch_lengths_;
  EigenVectorXd q_;
  EigenVectorXd unconditional_node_probabilities_;
  EigenVectorXd inverted_sbn_prior_;

  // The number of rows is equal to the number of GPCSPs.
  // The number of columns is equal to the number of site patterns.
  // The rows are indexed in the same way as branch_lengths_ and q_.
  // Entry (i,j) stores the marginal log likelihood over all trees that include
  // a GPCSP corresponding to index i at site j.
  EigenMatrixXd log_likelihoods_;

  // Internal "temporaries" useful for likelihood and derivative calculation.
  EigenVectorXd per_pattern_log_likelihoods_;
  EigenVectorXd per_pattern_likelihoods_;
  EigenVectorXd per_pattern_likelihood_derivatives_;
  EigenVectorXd per_pattern_likelihood_derivative_ratios_;
  EigenVectorXd per_pattern_likelihood_second_derivatives_;
  EigenVectorXd per_pattern_likelihood_second_derivative_ratios_;

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
  Eigen::Matrix4d hessian_matrix_;
  Eigen::Vector4d stationary_distribution_ = substitution_model_.GetFrequencies();
  EigenVectorXd site_pattern_weights_;

  // For hybrid marginal calculations. #328
  // The PLV coming down from the root.
  EigenMatrixXd quartet_root_plv_;
  // The R-PLV pointing leafward from s.
  EigenMatrixXd quartet_r_s_plv_;
  // The Q-PLV pointing leafward from s.
  EigenMatrixXd quartet_q_s_plv_;
  // The R-PLV pointing leafward from t.
  EigenMatrixXd quartet_r_sorted_plv_;

  void InitializePLVsWithSitePatterns();

  void RescalePLV(size_t plv_idx, int amount);
  void AssertPLVIsFinite(size_t plv_idx, const std::string& message) const;
  std::pair<double, double> PLVMinMax(size_t plv_idx) const;
  // If a PLV all entries smaller than rescaling_threshold_ then rescale it up and
  // increment the corresponding entry in rescaling_counts_.
  void RescalePLVIfNeeded(size_t plv_idx);
  double LogRescalingFor(size_t plv_idx);

  void TypeOfOptimization(const GPOperations::OptimizeBranchLength& op);
  void BrentOptimization(const GPOperations::OptimizeBranchLength& op);
  void GradientAscentOptimization(const GPOperations::OptimizeBranchLength& op);
  void LogSpaceGradientAscentOptimization(const GPOperations::OptimizeBranchLength& op);
  void NewtonOptimization(const GPOperations::OptimizeBranchLength& op);
  void AdaptiveGradientAscentOptimization(const GPOperations::OptimizeBranchLength& op);

  inline void PrepareUnrescaledPerPatternLikelihoodSecondDerivatives(size_t src1_idx,
                                                                     size_t src2_idx) {
    per_pattern_likelihood_second_derivatives_ =
        (plvs_.at(src1_idx).transpose() * hessian_matrix_ * plvs_.at(src2_idx))
            .diagonal()
            .array();
  }
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

  // This function is used to compute the marginal log likelihood over all trees that
  // have a given PCSP. We assume that transition_matrix_ is as desired, and src1_idx
  // and src2_idx are the two PLV indices on either side of the PCSP.
  inline void PreparePerPatternLogLikelihoodsForGPCSP(size_t src1_idx,
                                                      size_t src2_idx) {
    per_pattern_log_likelihoods_ =
        (plvs_.at(src1_idx).transpose() * transition_matrix_ * plvs_.at(src2_idx))
            .diagonal()
            .array()
            .log() +
        LogRescalingFor(src1_idx) + LogRescalingFor(src2_idx);
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("GPEngine") {
  EigenVectorXd empty_vector;
  SitePattern hello_site_pattern = SitePattern::HelloSitePattern();
  GPEngine engine(hello_site_pattern, 6 * 5, 5, "_ignore/mmapped_plv.data",
                  GPEngine::default_rescaling_threshold_, empty_vector, empty_vector,
                  empty_vector);
  engine.SetTransitionMatrixToHaveBranchLength(0.75);
  // Computed directly:
  // https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_%28Jukes_and_Cantor_1969%29
  CHECK(fabs(0.52590958087 - engine.GetTransitionMatrix()(0, 0)) < 1e-10);
  CHECK(fabs(0.1580301397 - engine.GetTransitionMatrix()(0, 1)) < 1e-10);
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_GP_ENGINE_HPP_
