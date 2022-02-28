// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A visitor for GPOperations. See
// https://arne-mertz.de/2018/05/modern-c-features-stdvariant-and-stdvisit/

#pragma once

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

  // Apply all operations in vector in order from beginning to end.
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

  size_t GetSitePatternCount() const { return site_pattern_.PatternCount(); };
  // Calculate a vector of likelihoods, one for each summand of the hybrid marginal.
  EigenVectorXd CalculateQuartetHybridLikelihoods(const QuartetHybridRequest& request);
  // Calculate the actual hybrid marginal and store it in the corresponding entry of
  // hybrid_marginal_log_likelihoods_.
  void ProcessQuartetHybridRequest(const QuartetHybridRequest& request);
  void PrintPLV(size_t plv_idx);
  // Use branch lengths from loaded sample as a starting point for optimization.
  void HotStartBranchLengths(const RootedTreeCollection& tree_collection,
                             const BitsetSizeMap& indexer);
  // Gather branch lengths from loaded sample with their corresponding pcsp.
  SizeDoubleVectorMap GatherBranchLengths(const RootedTreeCollection& tree_collection,
                                          const BitsetSizeMap& indexer);
  void FunctionOverRootedTreeCollection(
      std::function<void(size_t, const RootedTree&, const Node*)>
          function_on_tree_node_by_gpcsp,
      const RootedTreeCollection& tree_collection, const BitsetSizeMap& indexer);

  DoublePair LogLikelihoodAndDerivative(const GPOperations::OptimizeBranchLength& op);

  double PLVByteCount() const { return mmapped_master_plv_.ByteCount(); };

 private:
  //
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

 public:
  static constexpr double default_rescaling_threshold_ = 1e-40;
  // Initial branch length during first branch length opimization.
  static constexpr double default_branch_length_ = 0.1;

 private:
  // Absolute lower bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double min_log_branch_length_ = -13.9;
  // Absolute upper bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double max_log_branch_length_ = 1.1;
  // Precision used for checking convergence of branch length optimization.
  int significant_digits_for_optimization_ = 6;
  //
  double relative_tolerance_for_optimization_ = 1e-2;
  // Step size used for gradient-based branch length optimization.
  double step_size_for_optimization_ = 5e-4;
  // Number of iterations allowed for branch length optimization.
  size_t max_iter_for_optimization_ = 1000;

  // Descriptor containing all taxa and sequence alignments.
  SitePattern site_pattern_;
  // Rescaling threshold factor to prevent under/overflow errors.
  const double rescaling_threshold_;
  // Rescaling threshold in log space.
  const double log_rescaling_threshold_;

  // ** Per-Node Data

  // Total number of PLVs across entire DAG. Proportional to the number of nodes in DAG.
  // plv_count_ = node_count_without_root * plv_per_node.
  const size_t plv_count_per_node_ = 6;
  size_t plv_count_;
  // Master PLV: Large data block of virtual memory for Partial Likelihood Vectors.
  // Subdivided into sections for plvs_.
  std::unique_ptr<MmappedNucleotidePLV> mmapped_master_plv_ptr_ = nullptr;
  MmappedNucleotidePLV mmapped_master_plv_;
  // Partial Likelihood Vectors.
  // plvs_ store the following (see PLVHandler::GetPLVIndex):
  // [0, num_nodes): p(s).
  // [num_nodes, 2*num_nodes): phat(s_right).
  // [2*num_nodes, 3*num_nodes): phat(s_left).
  // [3*num_nodes, 4*num_nodes): rhat(s_right) = rhat(s_left).
  // [4*num_nodes, 5*num_nodes): r(s_right).
  // [5*num_nodes, 6*num_nodes): r(s_left).
  NucleotidePLVRefVector plvs_;
  // Rescaling count for each plv.
  EigenVectorXi rescaling_counts_;
  // For hybrid marginal calculations. #328
  // The PLV coming down from the root to s.
  EigenMatrixXd quartet_root_plv_;
  // The R-PLV pointing leafward from s.
  EigenMatrixXd quartet_r_s_plv_;
  // The Q-PLV pointing leafward from s.
  EigenMatrixXd quartet_q_s_plv_;
  // The R-PLV pointing leafward from t.
  EigenMatrixXd quartet_r_sorted_plv_;

  // ** Per-Edge Data

  // Total number of edges in DAG.
  size_t gpcsp_count_;
  // branch_lengths_, q_, etc. are indexed in the same way as sbn_parameters_ in
  // gp_instance.
  EigenVectorXd branch_lengths_;
  // During initialization, stores the SBN prior.
  // After UpdateSBNProbabilities(), stores the SBN probabilities.
  // Stored in log space.
  EigenVectorXd q_;
  EigenVectorXd unconditional_node_probabilities_;
  EigenVectorXd inverted_sbn_prior_;
  // The number of rows is equal to the number of GPCSPs.
  // The number of columns is equal to the number of site patterns.
  // The rows are indexed in the same way as branch_lengths_ and q_.
  // Entry (i,j) stores the marginal log likelihood over all trees that include
  // a GPCSP corresponding to index i at site j.
  EigenMatrixXd log_likelihoods_;
  // The length of this vector is equal to the number of site patterns.
  // Entry j stores the marginal log likelihood over all trees at site pattern
  // j.
  EigenVectorXd log_marginal_likelihood_;
  // This vector is indexed by the GPCSPs and stores the hybrid marginals if they are
  // available.
  EigenVectorXd hybrid_marginal_log_likelihoods_;
  // Internal "temporaries" useful for likelihood and derivative calculation.
  EigenVectorXd per_pattern_log_likelihoods_;
  EigenVectorXd per_pattern_likelihoods_;
  EigenVectorXd per_pattern_likelihood_derivatives_;
  EigenVectorXd per_pattern_likelihood_derivative_ratios_;

  // ** Model

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
