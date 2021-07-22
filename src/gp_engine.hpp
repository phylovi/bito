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
#include "reindexer.hpp"

class GPEngine {
 public:
  GPEngine(SitePattern site_pattern, size_t node_count, size_t plv_count_per_node,
           size_t gpcsp_count, const std::string& mmap_file_path,
           double rescaling_threshold, EigenVectorXd sbn_prior,
           EigenVectorXd unconditional_node_probabilities,
           EigenVectorXd inverted_sbn_prior);

  // Initialize prior with given starting values.
  void InitializePriors(EigenVectorXd sbn_prior,
                        EigenVectorXd unconditional_node_probabilities,
                        EigenVectorXd inverted_sbn_prior);
  // This sets all priors to 1.0. For testing purposes.
  void SetNullPrior();

  // ** Resizing and Reindexing

  // Resize GPEngine to accomodate DAG with given number of nodes and edges.  Option to
  // remap data according to DAG reindexers.  Option to give explicit numver of nodes or
  // edges to allocate memory for (this is the only way memory allocation will be
  // decreased).
  void GrowPLVs(const size_t node_count,
                std::optional<const Reindexer> node_reindexer = std::nullopt,
                std::optional<const size_t> explicit_allocation = std::nullopt,
                const bool on_initialization = false);
  void GrowGPCSPs(const size_t gpcsp_count,
                  std::optional<const Reindexer> gpcsp_reindexer = std::nullopt,
                  std::optional<const size_t> explicit_allocation = std::nullopt,
                  const bool on_intialization = false);
  // Remap node and edge-based data according to reordering of DAG nodes and edges.
  void ReindexPLVs(const Reindexer node_reindexer, const size_t old_node_count);
  void ReindexGPCSPs(const Reindexer gpcsp_reindexer, const size_t old_gpcsp_count);
  // Grow space for storing temporary computation.
  void GrowTempPLVs(const size_t new_node_padding);
  void GrowTempGPCSPs(const size_t new_gpcsp_padding);
  // Get index according to offset into temporary space.
  size_t GetTempPLVIndex(const size_t plv_offset) const;
  size_t GetTempGPCSPIndex(const size_t gpcsp_offset) const;

  // ** GPOperations

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
  void SetBranchLengths(EigenVectorXd branch_lengths);
  void SetBranchLengthsToConstant(double branch_length);
  void SetBranchLengthsToDefault();
  void ResetLogMarginalLikelihood();

  // The purpose of these functions are here to move data associated with the subsplit
  // DAG from their temporary locations before reindexing to their final locations after
  // reindexing (though they are more general).
  void CopyNodeData(const size_t src_node_idx, const size_t dest_node_idx);
  void CopyPLVData(const size_t src_plv_idx, const size_t dest_plv_idx);
  void CopyGPCSPData(const size_t src_gpcsp_idx, const size_t dest_gpcsp_idx);

  // ** Getters

  // Get Branch Lengths.
  EigenVectorXd GetBranchLengths() const;
  EigenVectorXd GetBranchLengths(const size_t start, const size_t length) const;
  // Get Branch Lengths from temporary space.
  EigenVectorXd GetTempBranchLengths(const size_t start, const size_t length) const;
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
  EigenVectorXd GetPerGPCSPLogLikelihoods(const size_t start,
                                          const size_t length = 1) const;
  // Get PerGPCSPLogLikelihoods from temporary space.
  EigenVectorXd GetTempPerGPCSPLogLikelihoods(const size_t start,
                                              const size_t length = 1) const;
  // This is the full marginal likelihood sum restricted to trees containing a PCSP.
  // When we sum the log of eq:PerGPCSPComponentsOfFullMarginal over the sites, we get
  // out a term that is the number of sites times the log of the prior conditional PCSP
  // probability.
  EigenVectorXd GetPerGPCSPComponentsOfFullLogMarginal() const;
  // #288 reconsider this name
  EigenConstMatrixXdRef GetLogLikelihoodMatrix() const;
  EigenConstVectorXdRef GetHybridMarginals() const;
  EigenConstVectorXdRef GetSBNParameters() const;

  double GetLogMarginalLikelihood() const;
  const Eigen::Matrix4d& GetTransitionMatrix() const { return transition_matrix_; };

  // Get Partial Likelihood Vector.
  EigenConstMatrixXdRef GetPLV(size_t plv_index) const { return plvs_.at(plv_index); }
  EigenConstMatrixXdRef GetTempPLV(size_t plv_index) const {
    return GetPLV(GetTempPLVIndex(plv_index));
  }

  // ** Other Operations

  // Calculate a vector of likelihoods, one for each summand of the hybrid marginal.
  EigenVectorXd CalculateQuartetHybridLikelihoods(const QuartetHybridRequest& request);
  // Calculate the actual hybrid marginal and store it in the corresponding entry of
  // hybrid_marginal_log_likelihoods_.
  void ProcessQuartetHybridRequest(const QuartetHybridRequest& request);

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

  // Take the first branch length encountered (in the supplied tree collection) for a
  // given edge for the branch length of the sDAG. Set branch lengths that are not thus
  // specified to default_branch_length_.
  void TakeFirstBranchLength(const RootedTreeCollection& tree_collection,
                             const BitsetSizeMap& indexer);

  DoublePair LogLikelihoodAndDerivative(const GPOperations::OptimizeBranchLength& op);
  std::tuple<double, double, double> LogLikelihoodAndFirstTwoDerivatives(
      const GPOperations::OptimizeBranchLength& op);

  // ** I/O

  // Output PLV to standard output.
  void PrintPLV(size_t plv_idx) const;
  // Output PLV to string.
  std::string PLVToString(size_t plv_idx) const;

  // ** Counts

  double PLVByteCount() const { return mmapped_master_plv_.ByteCount(); };
  size_t GetSitePatternCount() const { return site_pattern_.PatternCount(); };
  size_t GetNodeCount() const { return node_count_; };
  size_t GetTempNodeCount() const { return node_padding_; }
  size_t GetPaddedNodeCount() const { return node_count_ + node_padding_; };
  size_t GetAllocatedNodeCount() const { return node_alloc_; }
  size_t GetPLVCount() const { return GetNodeCount() * plv_count_per_node_; };
  size_t GetTempPLVCount() const { return GetTempNodeCount() * plv_count_per_node_; };
  size_t GetPaddedPLVCount() const {
    return GetPaddedNodeCount() * plv_count_per_node_;
  };
  size_t GetAllocatedPLVCount() const {
    return GetAllocatedNodeCount() * plv_count_per_node_;
  }
  size_t GetGPCSPCount() const { return gpcsp_count_; };
  size_t GetTempGPCSPCount() const { return gpcsp_padding_; };
  size_t GetPaddedGPCSPCount() const { return gpcsp_count_ + gpcsp_padding_; };
  size_t GetAllocatedGPCSPCount() const { return gpcsp_alloc_; };

 private:
  // Initialize PLVs and populate leaf PLVs with taxon site data.
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
  double relative_tolerance_for_optimization_ = 1e-3;
  double denominator_tolerance_for_newton_ = 1e-15;
  double step_size_for_optimization_ = 5e-4;
  double step_size_for_log_space_optimization_ = 1.0005;
  // Number of iterations allowed for branch length optimization.
  size_t max_iter_for_optimization_ = 1000;

  // Descriptor containing all taxa and sequence alignments.
  SitePattern site_pattern_;
  // Rescaling threshold factor to prevent under/overflow errors.
  const double rescaling_threshold_;
  // Rescaling threshold in log space.
  const double log_rescaling_threshold_;

  // ** Data Sizing
  // "Count" is the currently occupied by data.
  // "Padding" is the amount of free working space added to end of occupied space.
  // "Alloc" is the total current memory allocation.
  // "Resizing factor" is the amount of extra storage allocated for when resizing.

  // Total number of nodes in DAG.
  size_t node_count_ = 0;
  size_t node_alloc_ = 0;
  size_t node_padding_ = 2;
  // PLV count is proportional to the node_count.
  const size_t plv_count_per_node_ = 6;
  // Total number of edges in DAG. Determines sizes of data vectors stored on edges like
  // branch lengths.
  size_t gpcsp_count_ = 0;
  size_t gpcsp_alloc_ = 0;
  size_t gpcsp_padding_ = 3;
  // Growth factor when reallocating data.
  constexpr static double resizing_factor_ = 2.0;

  // ** Per-Node Data

  // Master PLV: Large data block of virtual memory for Partial Likelihood Vectors.
  // Subdivided into sections for plvs_.
  std::string mmap_file_path_;
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
  EigenVectorXd unconditional_node_probabilities_;
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

  // branch_lengths_, q_, etc. are indexed in the same way as sbn_parameters_ in
  // gp_instance.
  EigenVectorXd branch_lengths_;
  // During initialization, stores the SBN prior.
  // After UpdateSBNProbabilities(), stores the SBN probabilities.
  // Stored in log space.
  EigenVectorXd q_;
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
  EigenVectorXd per_pattern_likelihood_second_derivatives_;
  EigenVectorXd per_pattern_likelihood_second_derivative_ratios_;

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
  GPEngine engine(hello_site_pattern, 5, 6, 5, "_ignore/mmapped_plv.data",
                  GPEngine::default_rescaling_threshold_, empty_vector, empty_vector,
                  empty_vector);
  engine.SetTransitionMatrixToHaveBranchLength(0.75);
  // Computed directly:
  // https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_%28Jukes_and_Cantor_1969%29
  CHECK(fabs(0.52590958087 - engine.GetTransitionMatrix()(0, 0)) < 1e-10);
  CHECK(fabs(0.1580301397 - engine.GetTransitionMatrix()(0, 1)) < 1e-10);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
