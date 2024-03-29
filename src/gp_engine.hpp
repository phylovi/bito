// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A visitor for GPOperations. See
// https://arne-mertz.de/2018/05/modern-c-features-stdvariant-and-stdvisit/

#pragma once

#include "eigen_sugar.hpp"
#include "gp_operation.hpp"
#include "pv_handler.hpp"
#include "numerical_utils.hpp"
#include "quartet_hybrid_request.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "site_pattern.hpp"
#include "substitution_model.hpp"
#include "reindexer.hpp"
#include "subsplit_dag_storage.hpp"
#include "optimization.hpp"
#include "dag_branch_handler.hpp"
#include "dag_data.hpp"

class GPEngine {
 public:
  GPEngine(SitePattern site_pattern, size_t node_count, size_t gpcsp_count,
           const std::string& mmap_file_path, double rescaling_threshold,
           EigenVectorXd sbn_prior, EigenVectorXd unconditional_node_probabilities,
           EigenVectorXd inverted_sbn_prior, bool use_gradients);

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
  void ReindexPLVs(const Reindexer& node_reindexer, const size_t old_node_count);
  void ReindexGPCSPs(const Reindexer& gpcsp_reindexer, const size_t old_gpcsp_count);
  // Grow space for storing temporary computation.
  void GrowSparePLVs(const size_t new_node_spare_count);
  void GrowSpareGPCSPs(const size_t new_gpcsp_spare_count);

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

  // ** Branch Length Optimization

  void InitializeBranchLengthHandler();
  void OptimizeBranchLength(const GPOperations::OptimizeBranchLength& op);
  void SetOptimizationMethod(const OptimizationMethod method);
  void UseGradientOptimization(const bool use_gradients);
  void SetSignificantDigitsForOptimization(int significant_digits);
  size_t GetOptimizationCount() { return branch_handler_.GetOptimizationCount(); }
  void ResetOptimizationCount() { branch_handler_.ResetOptimizationCount(); }
  void IncrementOptimizationCount() { branch_handler_.IncrementOptimizationCount(); }
  bool IsFirstOptimization() { return branch_handler_.IsFirstOptimization(); }

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
  void CopyNodeData(const NodeId src_node_idx, const NodeId dest_node_idx);
  void CopyPLVData(const size_t src_plv_idx, const size_t dest_plv_idx);
  void CopyGPCSPData(const EdgeId src_gpcsp_idx, const EdgeId dest_gpcsp_idx);

  // ** Access

  // Get Branch Lengths.
  const DAGBranchHandler& GetBranchLengthHandler() const { return branch_handler_; }
  DAGBranchHandler& GetBranchLengthHandler() { return branch_handler_; }
  EigenVectorXd GetBranchLengths() const;
  EigenVectorXd GetBranchLengths(const size_t start, const size_t length) const;
  // Get Branch Lengths from temporary space.
  EigenVectorXd GetSpareBranchLengths(const size_t start, const size_t length) const;
  // Get differences for branch lengths during optimization to assess convergence.
  EigenVectorXd GetBranchLengthDifferences() const;

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
  EigenVectorXd GetSparePerGPCSPLogLikelihoods(const size_t start,
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

  // Partial Likelihood Vector Handler.
  const PLVNodeHandler& GetPLVHandler() const { return plv_handler_; }
  NucleotidePLVRefVector& GetPLVs() { return plv_handler_.GetPVs(); }
  const NucleotidePLVRefVector& GetPLVs() const { return plv_handler_.GetPVs(); }
  NucleotidePLVRef& GetPLV(const PVId plv_index) { return plv_handler_(plv_index); }
  const NucleotidePLVRef& GetPLV(const PVId plv_index) const {
    return plv_handler_(plv_index);
  }
  NucleotidePLVRef& GetSparePLV(const PVId plv_index) {
    return plv_handler_.GetSparePV(plv_index);
  }
  const NucleotidePLVRef& GetSparePLV(const PVId plv_index) const {
    return plv_handler_.GetSparePV(plv_index);
  }
  PVId GetSparePLVIndex(const PVId plv_index) const {
    return plv_handler_.GetSparePVIndex(plv_index);
  }

  // ** Other Operations

  // Calculate a vector of likelihoods, one for each summand of the hybrid marginal.
  EigenVectorXd CalculateQuartetHybridLikelihoods(const QuartetHybridRequest& request);
  // Calculate the actual hybrid marginal and store it in the corresponding entry of
  // hybrid_marginal_log_likelihoods_.
  void ProcessQuartetHybridRequest(const QuartetHybridRequest& request);

  // Gather branch lengths from loaded sample with their corresponding pcsp.
  SizeDoubleVectorMap GatherBranchLengths(const RootedTreeCollection& tree_collection,
                                          const BitsetSizeMap& indexer);

  // Use branch lengths from loaded sample as a starting point for optimization. Use the
  // mean branch length found for a given edge.
  void HotStartBranchLengths(const RootedTreeCollection& tree_collection,
                             const BitsetSizeMap& indexer);
  // Take the first branch length encountered (in the supplied tree collection) for a
  // given edge for the branch length of the sDAG. Set branch lengths that are not thus
  // specified to default_branch_length_.
  void TakeFirstBranchLength(const RootedTreeCollection& tree_collection,
                             const BitsetSizeMap& indexer);

  DoublePair LogLikelihoodAndDerivative(const GPOperations::OptimizeBranchLength& op);
  DoublePair LogLikelihoodAndDerivative(const size_t gpcsp, const size_t rootward,
                                        const size_t leafward);
  std::tuple<double, double, double> LogLikelihoodAndFirstTwoDerivatives(
      const GPOperations::OptimizeBranchLength& op);
  std::tuple<double, double, double> LogLikelihoodAndFirstTwoDerivatives(
      const size_t gpcsp, const size_t rootward, const size_t leafward);

  // ** I/O

  // Output PLV to string.
  std::string PLVToString(const PVId plv_idx) const;
  // Output LogLikelihood to string.
  std::string LogLikelihoodMatrixToString() const;

  // ** Counts

  size_t GetPLVCountPerNode() const { return plv_handler_.GetPVCountPer(); }
  size_t GetSitePatternCount() const { return site_pattern_.PatternCount(); };
  // Node Counts.
  size_t GetNodeCount() const { return plv_handler_.GetCount(); };
  size_t GetSpareNodeCount() const { return plv_handler_.GetSpareCount(); }
  size_t GetAllocatedNodeCount() const { return plv_handler_.GetAllocatedCount(); }
  size_t GetPaddedNodeCount() const { return plv_handler_.GetPaddedCount(); };
  void SetNodeCount(const size_t node_count) { plv_handler_.SetCount(node_count); }
  void SetSpareNodeCount(const size_t node_spare_count) {
    plv_handler_.SetSpareCount(node_spare_count);
  }
  void SetAllocatedNodeCount(const size_t node_alloc) {
    plv_handler_.SetAllocatedCount(node_alloc);
  }
  // PLV Counts.
  size_t GetPLVCount() const { return plv_handler_.GetPVCount(); };
  size_t GetSparePLVCount() const { return plv_handler_.GetSparePVCount(); };
  size_t GetPaddedPLVCount() const { return plv_handler_.GetPaddedPVCount(); };
  size_t GetAllocatedPLVCount() const { return plv_handler_.GetAllocatedPVCount(); }
  // GPCSP Counts.
  size_t GetGPCSPCount() const { return gpcsp_count_; };
  size_t GetSpareGPCSPCount() const { return gpcsp_spare_count_; };
  size_t GetAllocatedGPCSPCount() const { return gpcsp_alloc_; };
  size_t GetPaddedGPCSPCount() const { return GetGPCSPCount() + GetSpareGPCSPCount(); };
  size_t GetSpareGPCSPIndex(const size_t gpcsp_offset) const {
    const size_t gpcsp_scratch_size = GetPaddedGPCSPCount() - GetGPCSPCount();
    Assert(gpcsp_offset < gpcsp_scratch_size,
           "Requested gpcsp_offset outside of allocated scratch space.");
    return gpcsp_offset + GetGPCSPCount();
  }
  void SetGPCSPCount(const size_t gpcsp_count) { gpcsp_count_ = gpcsp_count; }
  void SetSpareGPCSPCount(const size_t gpcsp_spare_count) {
    gpcsp_spare_count_ = gpcsp_spare_count;
  }
  void SetAllocatedGPCSPCount(const size_t gpcsp_alloc) { gpcsp_alloc_ = gpcsp_alloc; }

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

  inline void PrepareUnrescaledPerPatternLikelihoodSecondDerivatives(size_t src1_idx,
                                                                     size_t src2_idx) {
    per_pattern_likelihood_second_derivatives_ =
        (GetPLV(PVId(src1_idx)).transpose() * hessian_matrix_ * GetPLV(PVId(src2_idx)))
            .diagonal()
            .array();
  }
  inline void PrepareUnrescaledPerPatternLikelihoodDerivatives(size_t src1_idx,
                                                               size_t src2_idx) {
    per_pattern_likelihood_derivatives_ = (GetPLV(PVId(src1_idx)).transpose() *
                                           derivative_matrix_ * GetPLV(PVId(src2_idx)))
                                              .diagonal()
                                              .array();
  }

  inline void PrepareUnrescaledPerPatternLikelihoods(size_t src1_idx, size_t src2_idx) {
    per_pattern_likelihoods_ = (GetPLV(PVId(src1_idx)).transpose() *
                                transition_matrix_ * GetPLV(PVId(src2_idx)))
                                   .diagonal()
                                   .array();
  }

  // This function is used to compute the marginal log likelihood over all trees that
  // have a given PCSP. We assume that transition_matrix_ is as desired, and src1_idx
  // and src2_idx are the two PLV indices on either side of the PCSP.
  inline void PreparePerPatternLogLikelihoodsForGPCSP(size_t src1_idx,
                                                      size_t src2_idx) {
    per_pattern_log_likelihoods_ = (GetPLV(PVId(src1_idx)).transpose() *
                                    transition_matrix_ * GetPLV(PVId(src2_idx)))
                                       .diagonal()
                                       .array()
                                       .log() +
                                   LogRescalingFor(src1_idx) +
                                   LogRescalingFor(src2_idx);
  }

 public:
  static constexpr double default_rescaling_threshold_ = 1e-40;

 private:
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
  // Note: All node and PLV counts are handled by the PLVNodeHandler.

  // Total number of edges in DAG. Determines sizes of data vectors indexed on edges
  // like branch lengths.
  size_t gpcsp_count_ = 0;
  size_t gpcsp_alloc_ = 0;
  size_t gpcsp_spare_count_ = 3;
  // Growth factor when reallocating data.
  constexpr static double resizing_factor_ = 2.0;

  // ** Per-Node Data

  // Partial Likelihood Vector Handler.
  PLVNodeHandler plv_handler_;
  // Unconditional probabilites for each node in DAG.
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

  // branch_handler_, q_, etc. are indexed in the same way as sbn_parameters_ in
  // gp_instance.
  DAGBranchHandler branch_handler_;
  // During initialization, stores the SBN prior.
  // After UpdateSBNProbabilities(), stores the SBN probabilities.
  // Stored in log space.
  EigenVectorXd q_;
  EigenVectorXd inverted_sbn_prior_;

  // The number of rows is equal to the number of GPCSPs.
  // The number of columns is equal to the number of site patterns.
  // The rows are indexed in the same way as branch_handler_ and q_.
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
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("GPEngine") {
  EigenVectorXd empty_vector;
  SitePattern hello_site_pattern = SitePattern::HelloSitePattern();
  GPEngine engine(hello_site_pattern, 5, 5, "_ignore/mmapped_plv.data",
                  GPEngine::default_rescaling_threshold_, empty_vector, empty_vector,
                  empty_vector, false);
  engine.SetTransitionMatrixToHaveBranchLength(0.75);
  // Computed directly:
  // https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_%28Jukes_and_Cantor_1969%29
  CHECK(fabs(0.52590958087 - engine.GetTransitionMatrix()(0, 0)) < 1e-10);
  CHECK(fabs(0.1580301397 - engine.GetTransitionMatrix()(0, 1)) < 1e-10);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
