// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"

#include "optimization.hpp"
#include "sugar.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t node_count,
                   size_t plv_count_per_node, size_t gpcsp_count,
                   const std::string& mmap_file_path, double rescaling_threshold,
                   EigenVectorXd sbn_prior,
                   EigenVectorXd unconditional_node_probabilities,
                   EigenVectorXd inverted_sbn_prior, bool use_gradients)
    : site_pattern_(std::move(site_pattern)),
      rescaling_threshold_(rescaling_threshold),
      log_rescaling_threshold_(log(rescaling_threshold)),
      plv_handler_(mmap_file_path, 0, site_pattern_.PatternCount(), resizing_factor_),
      unconditional_node_probabilities_(std::move(unconditional_node_probabilities)),
      q_(std::move(sbn_prior)),
      inverted_sbn_prior_(std::move(inverted_sbn_prior)) {
  // Initialize site pattern-based data
  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ = EigenVectorXdOfStdVectorDouble(weights);
  log_marginal_likelihood_.resize(site_pattern_.PatternCount());
  log_marginal_likelihood_.setConstant(DOUBLE_NEG_INF);
  // Initialize node-based data
  GrowPLVs(node_count, std::nullopt, std::nullopt, true);
  InitializePLVsWithSitePatterns();
  // Initialize edge-based data
  GrowGPCSPs(gpcsp_count, std::nullopt, std::nullopt, true);
  // Initialize PLV temporaries.
  quartet_root_plv_ = plv_handler_.GetPLV(0);
  quartet_root_plv_.setZero();
  quartet_r_s_plv_ = quartet_root_plv_;
  quartet_q_s_plv_ = quartet_root_plv_;
  quartet_r_sorted_plv_ = quartet_root_plv_;

  optimization_method_ =
      (use_gradients ? OptimizationMethod::BrentOptimizationWithGradients
                     : OptimizationMethod::BrentOptimization);
}

void GPEngine::InitializePriors(EigenVectorXd sbn_prior,
                                EigenVectorXd unconditional_node_probabilities,
                                EigenVectorXd inverted_sbn_prior) {
  Assert(size_t(unconditional_node_probabilities.size()) == GetNodeCount(),
         "unconditional_node_probabilities is wrong size for GPEngine.");
  Assert(size_t(sbn_prior.size()) == GetGPCSPCount(),
         "sbn_prior is wrong size for GPEngine.");
  Assert(size_t(inverted_sbn_prior.size()) == GetGPCSPCount(),
         "inverted_sbn_prior is wrong size for GPEngine.");
  unconditional_node_probabilities_.segment(0, GetNodeCount()) =
      unconditional_node_probabilities;
  q_.segment(0, GetGPCSPCount()) = sbn_prior;
  inverted_sbn_prior_.segment(0, GetGPCSPCount()) = inverted_sbn_prior;
}

void GPEngine::SetNullPrior() { q_.setConstant(1.0); }

// ** Resize and Reindex

void GPEngine::GrowPLVs(const size_t new_node_count,
                        std::optional<const Reindexer> node_reindexer,
                        std::optional<const size_t> explicit_allocation,
                        const bool on_initialization) {
  const size_t old_node_count = GetNodeCount();
  const size_t old_plv_count = GetPLVCount();
  SetNodeCount(new_node_count);
  // Reallocate more space if needed.
  if ((GetPaddedNodeCount() > GetAllocatedNodeCount()) ||
      explicit_allocation.has_value()) {
    SetAllocatedNodeCount(
        size_t(ceil(double(GetPaddedNodeCount()) * resizing_factor_)));
    if (explicit_allocation.has_value()) {
      Assert(explicit_allocation.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      SetAllocatedNodeCount(explicit_allocation.value() + GetTempNodeCount());
    }
    plv_handler_.Resize(new_node_count, GetAllocatedNodeCount());
    rescaling_counts_.conservativeResize(GetAllocatedPLVCount());
    unconditional_node_probabilities_.conservativeResize(GetAllocatedNodeCount());
  }
  // Resize to fit without deallocating unused memory.
  rescaling_counts_.conservativeResize(GetPaddedPLVCount());
  if (on_initialization) {
    unconditional_node_probabilities_.conservativeResize(GetNodeCount());
  } else {
    unconditional_node_probabilities_.conservativeResize(GetPaddedNodeCount());
  }
  // Initialize new work space.
  Assert((plv_handler_.GetPLVs().back().rows() == MmappedNucleotidePLV::base_count_) &&
             (plv_handler_.GetPLVs().back().cols() ==
              static_cast<Eigen::Index>(site_pattern_.PatternCount())) &&
             (size_t(plv_handler_.GetPLVs().size()) == GetAllocatedPLVCount()),
         "Didn't get the right shape of PLVs out of Subdivide.");
  for (size_t i = old_plv_count; i < GetPaddedPLVCount(); i++) {
    rescaling_counts_[i] = 0;
    plv_handler_.GetPLV(i).setZero();
  }
  for (size_t i = old_node_count; i < GetNodeCount(); i++) {
    if (!on_initialization) {
      unconditional_node_probabilities_[i] = 1.;
    }
  }
  // Reindex work space to realign with DAG.
  if (node_reindexer.has_value()) {
    ReindexPLVs(node_reindexer.value(), old_node_count);
  }
}

void GPEngine::GrowGPCSPs(const size_t new_gpcsp_count,
                          std::optional<const Reindexer> gpcsp_reindexer,
                          std::optional<const size_t> explicit_allocation,
                          const bool on_initialization) {
  const size_t old_gpcsp_count = GetGPCSPCount();
  SetGPCSPCount(new_gpcsp_count);
  // Reallocate more space if needed.
  if ((GetPaddedGPCSPCount() > GetAllocatedGPCSPCount()) ||
      explicit_allocation.has_value()) {
    SetAllocatedGPCSPCount(
        size_t(ceil(double(GetPaddedGPCSPCount()) * resizing_factor_)));
    if (explicit_allocation.has_value()) {
      Assert(explicit_allocation.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      SetAllocatedGPCSPCount(explicit_allocation.value() + GetTempGPCSPCount());
    }
    branch_lengths_.conservativeResize(GetAllocatedGPCSPCount());
    branch_length_differences_.conservativeResize(GetAllocatedGPCSPCount());
    hybrid_marginal_log_likelihoods_.conservativeResize(GetAllocatedGPCSPCount());
    log_likelihoods_.conservativeResize(GetAllocatedGPCSPCount(),
                                        site_pattern_.PatternCount());
    q_.conservativeResize(GetAllocatedGPCSPCount());
    inverted_sbn_prior_.conservativeResize(GetAllocatedGPCSPCount());
  }
  // Resize to fit without deallocating unused memory.
  branch_lengths_.conservativeResize(GetPaddedGPCSPCount());
  branch_length_differences_.conservativeResize(GetPaddedGPCSPCount());
  hybrid_marginal_log_likelihoods_.conservativeResize(GetPaddedGPCSPCount());
  log_likelihoods_.conservativeResize(GetPaddedGPCSPCount(),
                                      site_pattern_.PatternCount());

  if (on_initialization) {
    q_.conservativeResize(GetGPCSPCount());
    inverted_sbn_prior_.conservativeResize(GetGPCSPCount());
  } else {
    q_.conservativeResize(GetPaddedGPCSPCount());
    inverted_sbn_prior_.conservativeResize(GetPaddedGPCSPCount());
  }
  // Initialize new work space.
  for (size_t i = old_gpcsp_count; i < GetPaddedGPCSPCount(); i++) {
    branch_lengths_[i] = default_branch_length_;
    branch_length_differences_[i] = default_branch_length_;
    hybrid_marginal_log_likelihoods_[i] = DOUBLE_NEG_INF;
  }
  if (on_initialization) {
    branch_length_differences_[old_gpcsp_count] = 0;
  } else {
    for (size_t i = old_gpcsp_count; i < GetGPCSPCount(); i++) {
      q_[i] = 1.;
      inverted_sbn_prior_[i] = 1.;
    }
  }
  // Reindex work space to realign with DAG.
  if (gpcsp_reindexer.has_value()) {
    ReindexGPCSPs(gpcsp_reindexer.value(), old_gpcsp_count);
  }
}

void GPEngine::ReindexPLVs(const Reindexer& node_reindexer,
                           const size_t old_node_count) {
  Assert(node_reindexer.size() == GetNodeCount(),
         "Node Reindexer is the wrong size for GPEngine.");
  Assert(node_reindexer.IsValid(GetNodeCount()), "Node Reindexer is not valid.");

  // Expand node_reindexer into plv_reindexer.
  Reindexer plv_reindexer =
      plv_handler_.BuildPLVReindexer(node_reindexer, old_node_count, GetNodeCount());
  // Reindex data vectors
  plv_handler_.Reindex(plv_reindexer);
  Reindexer::ReindexInPlace<EigenVectorXi, int>(rescaling_counts_, plv_reindexer,
                                                plv_handler_.GetPLVCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(unconditional_node_probabilities_,
                                                   node_reindexer, GetNodeCount());
}

void GPEngine::ReindexGPCSPs(const Reindexer& gpcsp_reindexer,
                             const size_t old_gpcsp_count) {
  Assert(gpcsp_reindexer.size() == GetGPCSPCount(),
         "GPCSP Reindexer is the wrong size for GPEngine.");
  Assert(gpcsp_reindexer.IsValid(GetGPCSPCount()),
         "GPCSP Reindexer is not valid for GPEngine size.");
  // Reindex data vectors.
  Reindexer::ReindexInPlace<EigenVectorXd, double>(branch_lengths_, gpcsp_reindexer,
                                                   GetGPCSPCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(branch_length_differences_,
                                                   gpcsp_reindexer, GetGPCSPCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(hybrid_marginal_log_likelihoods_,
                                                   gpcsp_reindexer, GetGPCSPCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(q_, gpcsp_reindexer,
                                                   GetGPCSPCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(inverted_sbn_prior_, gpcsp_reindexer,
                                                   GetGPCSPCount());
}

void GPEngine::GrowTempPLVs(const size_t new_node_padding) {
  if (new_node_padding > GetTempNodeCount()) {
    SetTempNodeCount(new_node_padding);
    GrowPLVs(GetNodeCount());
  }
}

void GPEngine::GrowTempGPCSPs(const size_t new_gpcsp_padding) {
  if (new_gpcsp_padding > GetTempGPCSPCount()) {
    SetTempGPCSPCount(new_gpcsp_padding);
    GrowGPCSPs(GetGPCSPCount());
  }
}

size_t GPEngine::GetTempPLVIndex(const size_t plv_offset) const {
  const size_t plv_scratch_size = GetPaddedPLVCount() - GetPLVCount();
  Assert(plv_offset < plv_scratch_size,
         "Requested plv_offset outside of allocated scratch space.");
  return plv_offset + GetPLVCount();
}

size_t GPEngine::GetTempGPCSPIndex(const size_t gpcsp_offset) const {
  const size_t gpcsp_scratch_size = GetPaddedGPCSPCount() - GetGPCSPCount();
  Assert(gpcsp_offset < gpcsp_scratch_size,
         "Requested gpcsp_offset outside of allocated scratch space.");
  return gpcsp_offset + GetGPCSPCount();
}

// ** GPOperations

void GPEngine::operator()(const GPOperations::ZeroPLV& op) {
  plv_handler_.GetPLV(op.dest_).setZero();
  rescaling_counts_(op.dest_) = 0;
}

void GPEngine::operator()(const GPOperations::SetToStationaryDistribution& op) {
  auto& plv = plv_handler_.GetPLV(op.dest_);
  for (Eigen::Index row_idx = 0; row_idx < plv.rows(); ++row_idx) {
    // Multiplication by q_ avoids special treatment of the rhat vector for the
    // rootsplits.
    plv.row(row_idx).array() =
        q_(op.root_gpcsp_idx_) * stationary_distribution_(row_idx);
  }
  rescaling_counts_(op.dest_) = 0;
}

void GPEngine::operator()(const GPOperations::IncrementWithWeightedEvolvedPLV& op) {
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(op.gpcsp_));
  // We assume that we've done a PrepForMarginalization operation, and thus the
  // rescaling count for op.dest_ is the minimum of the rescaling counts among the
  // op.src_s. Thus this should be non-negative:
  const int rescaling_difference =
      rescaling_counts_(op.src_) - rescaling_counts_(op.dest_);
  Assert(rescaling_difference >= 0,
         "dest_ rescaling too large in IncrementWithWeightedEvolvedPLV");
  const double rescaling_factor =
      rescaling_difference == 0
          ? 1.
          : pow(rescaling_threshold_, static_cast<double>(rescaling_difference));
  // We are going to have evidence of reduced-precision arithmetic here because we are
  // adding together things of radically different rescaling amounts. This appears
  // unavoidable without special-purpose truncation code, which doesn't seem
  // worthwhile.
  plv_handler_.GetPLV(op.dest_) += rescaling_factor * q_(op.gpcsp_) *
                                   transition_matrix_ * plv_handler_.GetPLV(op.src_);
}

void GPEngine::operator()(const GPOperations::ResetMarginalLikelihood& op) {  // NOLINT
  ResetLogMarginalLikelihood();
}

void GPEngine::operator()(const GPOperations::IncrementMarginalLikelihood& op) {
  Assert(rescaling_counts_(op.stationary_times_prior_) == 0,
         "Surprise! Rescaled stationary distribution in IncrementMarginalLikelihood");
  // This operation does two things: imcrement the overall per-site log marginal
  // likelihood, and also set the conditional per-rootsplit marginal likelihood.
  //
  // We first calculate the unconditional contribution of the rootsplit to the overall
  // per-site marginal likelihood. It's an unconditional contribution because our
  // stationary distribution incorporates the prior on rootsplits.
  log_likelihoods_.row(op.rootsplit_) =
      (plv_handler_.GetPLV(op.stationary_times_prior_).transpose() *
       plv_handler_.GetPLV(op.p_))
          .diagonal()
          .array()
          .log() +
      LogRescalingFor(op.p_);
  // We can then increment the overall per-site marginal likelihood.
  log_marginal_likelihood_ = NumericalUtils::LogAddVectors(
      log_marginal_likelihood_, log_likelihoods_.row(op.rootsplit_));
  // However, we want the row in log_likelihoods_ to be the marginal likelihood
  // *conditional* on that rootsplit, so we log-divide by the rootsplit's probability.
  log_likelihoods_.row(op.rootsplit_).array() -= log(q_[op.rootsplit_]);
}

void GPEngine::operator()(const GPOperations::Multiply& op) {
  plv_handler_.GetPLV(op.dest_).array() =
      plv_handler_.GetPLV(op.src1_).array() * plv_handler_.GetPLV(op.src2_).array();
  rescaling_counts_(op.dest_) =
      rescaling_counts_(op.src1_) + rescaling_counts_(op.src2_);
  AssertPLVIsFinite(op.dest_, "Multiply dest_ is not finite");
  RescalePLVIfNeeded(op.dest_);
}

void GPEngine::operator()(const GPOperations::Likelihood& op) {
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(op.dest_));
  PreparePerPatternLogLikelihoodsForGPCSP(op.parent_, op.child_);
  log_likelihoods_.row(op.dest_) = per_pattern_log_likelihoods_;
}

void GPEngine::operator()(const GPOperations::OptimizeBranchLength& op) {
  return Optimization(op);
}

EigenVectorXd NormalizedPosteriorOfLogUnnormalized(
    EigenVectorXd log_unnormalized_posterior) {
  const double log_norm = NumericalUtils::LogSum(log_unnormalized_posterior);
  log_unnormalized_posterior.array() -= log_norm;
  return log_unnormalized_posterior.array().exp();
}

void GPEngine::operator()(const GPOperations::UpdateSBNProbabilities& op) {
  const size_t range_length = op.stop_ - op.start_;
  if (range_length == 1) {
    q_(op.start_) = 1.;
  } else {
    EigenVectorXd log_likelihoods;
    EigenConstVectorXdRef our_hybrid_log_likelihoods =
        hybrid_marginal_log_likelihoods_.segment(op.start_, range_length);
    if (our_hybrid_log_likelihoods.minCoeff() > DOUBLE_NEG_INF) {
      log_likelihoods = our_hybrid_log_likelihoods;
    } else {
      log_likelihoods = GetPerGPCSPLogLikelihoods(op.start_, range_length);
    }
    EigenVectorXd log_prior = q_.segment(op.start_, range_length).array().log();
    q_.segment(op.start_, range_length) =
        NormalizedPosteriorOfLogUnnormalized(log_likelihoods + log_prior);
  }
}

void GPEngine::operator()(const GPOperations::PrepForMarginalization& op) {
  const size_t src_count = op.src_vector_.size();
  Assert(src_count > 0, "Empty src_vector in PrepForMarginalization");
  SizeVector src_rescaling_counts(src_count);
  for (size_t idx = 0; idx < src_count; ++idx) {
    src_rescaling_counts[idx] = rescaling_counts_(op.src_vector_[idx]);
  }
  const auto min_rescaling_count =
      *std::min_element(src_rescaling_counts.begin(), src_rescaling_counts.end());
  rescaling_counts_(op.dest_) = min_rescaling_count;
}

void GPEngine::ProcessOperations(GPOperationVector operations) {
  for (const auto& operation : operations) {
    std::visit(*this, operation);
  }
}

void GPEngine::SetTransitionMatrixToHaveBranchLength(double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void GPEngine::SetTransitionAndDerivativeMatricesToHaveBranchLength(
    double branch_length) {
  diagonal_vector_ = (branch_length * eigenvalues_).array().exp();
  diagonal_matrix_.diagonal() = diagonal_vector_;
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
  // Now calculating derivative matrix
  diagonal_matrix_.diagonal() = eigenvalues_.array() * diagonal_vector_.array();
  derivative_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
  // Now calculating hessian matrix
  diagonal_matrix_.diagonal() =
      eigenvalues_.array() * eigenvalues_.array() * diagonal_vector_.array();
  hessian_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void GPEngine::SetTransitionMatrixToHaveBranchLengthAndTranspose(double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ =
      inverse_eigenmatrix_.transpose() * diagonal_matrix_ * eigenmatrix_.transpose();
}

void GPEngine::SetBranchLengths(EigenVectorXd branch_lengths) {
  Assert(size_t(branch_lengths.size()) == GetGPCSPCount(),
         "Size mismatch in GPEngine::SetBranchLengths.");
  branch_lengths_.segment(0, GetGPCSPCount()) = branch_lengths;
};

void GPEngine::SetBranchLengthsToConstant(double branch_length) {
  branch_lengths_.setConstant(branch_length);
};

void GPEngine::SetBranchLengthsToDefault() {
  branch_lengths_.setConstant(default_branch_length_);
};

void GPEngine::ResetLogMarginalLikelihood() {
  log_marginal_likelihood_.setConstant(DOUBLE_NEG_INF);
}

void GPEngine::CopyNodeData(const size_t src_node_idx, const size_t dest_node_idx) {
  Assert(
      (src_node_idx < GetPaddedNodeCount()) && (dest_node_idx < GetPaddedNodeCount()),
      "Cannot copy node data with src or dest index out-of-range.");
  unconditional_node_probabilities_[dest_node_idx] =
      unconditional_node_probabilities_[src_node_idx];
}

void GPEngine::CopyPLVData(const size_t src_plv_idx, const size_t dest_plv_idx) {
  Assert((src_plv_idx < GetPaddedPLVCount()) && (dest_plv_idx < GetPaddedPLVCount()),
         "Cannot copy PLV data with src or dest index out-of-range.");
  plv_handler_.GetPLV(dest_plv_idx) = plv_handler_.GetPLV(src_plv_idx);
  rescaling_counts_[dest_plv_idx] = rescaling_counts_[src_plv_idx];
}

void GPEngine::CopyGPCSPData(const size_t src_gpcsp_idx, const size_t dest_gpcsp_idx) {
  Assert((src_gpcsp_idx < GetPaddedGPCSPCount()) &&
             (dest_gpcsp_idx < GetPaddedGPCSPCount()),
         "Cannot copy PLV data with src or dest index out-of-range.");
  branch_lengths_[dest_gpcsp_idx] = branch_lengths_[src_gpcsp_idx];
  q_[dest_gpcsp_idx] = branch_lengths_[src_gpcsp_idx];
  inverted_sbn_prior_[dest_gpcsp_idx] = inverted_sbn_prior_[src_gpcsp_idx];
}

// ** Getters

double GPEngine::GetLogMarginalLikelihood() const {
  return (log_marginal_likelihood_.array() * site_pattern_weights_.array()).sum();
}

EigenVectorXd GPEngine::GetBranchLengths() const {
  return branch_lengths_.segment(0, GetGPCSPCount());
};

EigenVectorXd GPEngine::GetBranchLengths(const size_t start,
                                         const size_t length) const {
  Assert(start + length <= GetPaddedGPCSPCount(),
         "Requested range of BranchLengths is out-of-range.");
  return branch_lengths_.segment(start, length);
};

EigenVectorXd GPEngine::GetTempBranchLengths(const size_t start,
                                             const size_t length) const {
  return GetBranchLengths(GetTempGPCSPIndex(start), length);
}

EigenVectorXd GPEngine::GetBranchLengthDifferences() const {
  return branch_length_differences_.segment(0, GetGPCSPCount());
};

EigenVectorXd GPEngine::GetPerGPCSPLogLikelihoods() const {
  return log_likelihoods_.block(0, 0, GetGPCSPCount(), log_likelihoods_.cols()) *
         site_pattern_weights_;
};

EigenVectorXd GPEngine::GetPerGPCSPLogLikelihoods(const size_t start,
                                                  const size_t length) const {
  Assert(start + length <= GetPaddedGPCSPCount(),
         "Requested range of PerGPCSPLogLikelihoods is out-of-range.");
  return log_likelihoods_.block(start, 0, length, log_likelihoods_.cols()) *
         site_pattern_weights_;
};

EigenVectorXd GPEngine::GetTempPerGPCSPLogLikelihoods(const size_t start,
                                                      const size_t length) const {
  return GetPerGPCSPLogLikelihoods(GetTempGPCSPIndex(start), length);
}

EigenVectorXd GPEngine::GetPerGPCSPComponentsOfFullLogMarginal() const {
  return GetPerGPCSPLogLikelihoods().array() +
         static_cast<double>(site_pattern_.SiteCount()) * q_.array().log();
};

EigenConstMatrixXdRef GPEngine::GetLogLikelihoodMatrix() const {
  return log_likelihoods_.block(0, 0, GetGPCSPCount(), log_likelihoods_.cols());
};

EigenConstVectorXdRef GPEngine::GetHybridMarginals() const {
  return hybrid_marginal_log_likelihoods_;
};

EigenConstVectorXdRef GPEngine::GetSBNParameters() const { return q_; };

// ** Other Operations

DoublePair GPEngine::LogLikelihoodAndDerivative(
    const GPOperations::OptimizeBranchLength& op) {
  SetTransitionAndDerivativeMatricesToHaveBranchLength(branch_lengths_(op.gpcsp_));
  PreparePerPatternLogLikelihoodsForGPCSP(op.rootward_, op.leafward_);
  // The prior is expressed using the current value of q_.
  // The phylogenetic component of the likelihood is weighted with the number of times
  // we see the site patterns.
  const double log_likelihood = per_pattern_log_likelihoods_.dot(site_pattern_weights_);

  // The per-site likelihood derivative is calculated in the same way as the per-site
  // likelihood, but using the derivative matrix instead of the transition matrix.
  // We first prepare two useful vectors _without_ likelihood rescaling, because the
  // rescalings cancel out in the ratio below.
  PrepareUnrescaledPerPatternLikelihoodDerivatives(op.rootward_, op.leafward_);
  PrepareUnrescaledPerPatternLikelihoods(op.rootward_, op.leafward_);
  // If l_i is the per-site likelihood, the derivative of log(l_i) is the derivative
  // of l_i divided by l_i.
  per_pattern_likelihood_derivative_ratios_ =
      per_pattern_likelihood_derivatives_.array() / per_pattern_likelihoods_.array();
  const double log_likelihood_derivative =
      per_pattern_likelihood_derivative_ratios_.dot(site_pattern_weights_);
  return {log_likelihood, log_likelihood_derivative};
}

std::tuple<double, double, double> GPEngine::LogLikelihoodAndFirstTwoDerivatives(
    const GPOperations::OptimizeBranchLength& op) {
  SetTransitionAndDerivativeMatricesToHaveBranchLength(branch_lengths_(op.gpcsp_));
  PreparePerPatternLogLikelihoodsForGPCSP(op.rootward_, op.leafward_);

  const double log_likelihood = per_pattern_log_likelihoods_.dot(site_pattern_weights_);

  // The per-site likelihood derivative is calculated in the same way as the per-site
  // likelihood, but using the derivative matrix instead of the transition matrix.
  // We first prepare two useful vectors _without_ likelihood rescaling, because the
  // rescalings cancel out in the ratio below.
  PrepareUnrescaledPerPatternLikelihoodDerivatives(op.rootward_, op.leafward_);
  PrepareUnrescaledPerPatternLikelihoods(op.rootward_, op.leafward_);
  // If l_i is the per-site likelihood, the derivative of log(l_i) is the derivative
  // of l_i divided by l_i.
  per_pattern_likelihood_derivative_ratios_ =
      per_pattern_likelihood_derivatives_.array() / per_pattern_likelihoods_.array();
  const double log_likelihood_gradient =
      per_pattern_likelihood_derivative_ratios_.dot(site_pattern_weights_);

  // Second derivative is calculated the same way, but has an extra term due to
  // the product rule.

  PrepareUnrescaledPerPatternLikelihoodSecondDerivatives(op.rootward_, op.leafward_);

  per_pattern_likelihood_second_derivative_ratios_ =
      (per_pattern_likelihood_second_derivatives_.array() *
           per_pattern_likelihoods_.array() -
       per_pattern_likelihood_derivatives_.array() *
           per_pattern_likelihood_derivatives_.array()) /
      (per_pattern_likelihoods_.array() * per_pattern_likelihoods_.array());

  const double log_likelihood_hessian =
      per_pattern_likelihood_second_derivative_ratios_.dot(site_pattern_weights_);

  return std::make_tuple(log_likelihood, log_likelihood_gradient,
                         log_likelihood_hessian);
}

void GPEngine::InitializePLVsWithSitePatterns() {
  for (auto& plv : plv_handler_.GetPLVs()) {
    plv.setZero();
  }
  size_t taxon_idx = 0;
  for (const auto& pattern : site_pattern_.GetPatterns()) {
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative symbol!");
      if (symbol == MmappedNucleotidePLV::base_count_) {  // Gap character.
        plv_handler_.GetPLV(taxon_idx).col(site_idx).setConstant(1.);
      } else if (symbol < MmappedNucleotidePLV::base_count_) {
        plv_handler_.GetPLV(taxon_idx)(symbol, site_idx) = 1.;
      }
      site_idx++;
    }
    taxon_idx++;
  }
}

void GPEngine::RescalePLV(size_t plv_idx, int rescaling_count) {
  if (rescaling_count == 0) {
    return;
  }
  // else
  Assert(rescaling_count >= 0, "Negative rescaling count in RescalePLV.");
  plv_handler_.GetPLV(plv_idx) /=
      pow(rescaling_threshold_, static_cast<double>(rescaling_count));
  rescaling_counts_(plv_idx) += rescaling_count;
}

void GPEngine::AssertPLVIsFinite(size_t plv_idx, const std::string& message) const {
  Assert(plv_handler_.GetPLV(plv_idx).array().isFinite().all(), message);
}

std::pair<double, double> GPEngine::PLVMinMax(size_t plv_idx) const {
  return {plv_handler_.GetPLV(plv_idx).minCoeff(),
          plv_handler_.GetPLV(plv_idx).maxCoeff()};
}

void GPEngine::RescalePLVIfNeeded(size_t plv_idx) {
  auto [min_entry, max_entry] = PLVMinMax(plv_idx);
  Assert(min_entry >= 0., "PLV with negative entry (" + std::to_string(min_entry) +
                              ") passed to RescalePLVIfNeeded");
  if (max_entry == 0) {
    return;
  }
  // else
  int rescaling_count = 0;
  while (max_entry < rescaling_threshold_) {
    max_entry /= rescaling_threshold_;
    rescaling_count++;
  }
  RescalePLV(plv_idx, rescaling_count);
}

double GPEngine::LogRescalingFor(size_t plv_idx) {
  return static_cast<double>(rescaling_counts_(plv_idx)) * log_rescaling_threshold_;
}

void GPEngine::SetOptimizationMethod(const GPEngine::OptimizationMethod method) {
  optimization_method_ = method;
}

void GPEngine::Optimization(const GPOperations::OptimizeBranchLength& op) {
  switch (optimization_method_) {
    case OptimizationMethod::BrentOptimization:
      return BrentOptimization(op);
    case OptimizationMethod::BrentOptimizationWithGradients:
      return BrentOptimizationWithGradients(op);
    case OptimizationMethod::GradientAscentOptimization:
      return GradientAscentOptimization(op);
    case OptimizationMethod::LogSpaceGradientAscentOptimization:
      return LogSpaceGradientAscentOptimization(op);
    case OptimizationMethod::NewtonOptimization:
      return NewtonOptimization(op);
    default:
      Failwith("GPEngine::Optimization(): Invalid OptimizationMethod given.");
  }
}

void GPEngine::SetSignificantDigitsForOptimization(int significant_digits) {
  significant_digits_for_optimization_ = significant_digits;
}

void GPEngine::BrentOptimization(const GPOperations::OptimizeBranchLength& op) {
  if (branch_length_differences_(op.gpcsp_) < branch_length_difference_threshold_) {
    return;
  }

  auto negative_log_likelihood = [this, &op](double log_branch_length) {
    SetTransitionMatrixToHaveBranchLength(exp(log_branch_length));
    PreparePerPatternLogLikelihoodsForGPCSP(op.rootward_, op.leafward_);
    return -per_pattern_log_likelihoods_.dot(site_pattern_weights_);
  };

  double current_log_branch_length = log(branch_lengths_(op.gpcsp_));
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
    branch_lengths_(op.gpcsp_) = exp(current_log_branch_length);
  } else {
    branch_lengths_(op.gpcsp_) = exp(log_branch_length);
  }
  branch_length_differences_(op.gpcsp_) =
      abs(exp(current_log_branch_length) - branch_lengths_(op.gpcsp_));
}

void GPEngine::BrentOptimizationWithGradients(
    const GPOperations::OptimizeBranchLength& op) {
  if (branch_length_differences_(op.gpcsp_) < branch_length_difference_threshold_) {
    return;
  }

  auto negative_log_likelihood_and_derivative = [this, &op](double log_branch_length) {
    double branch_length = exp(log_branch_length);
    branch_lengths_(op.gpcsp_) = branch_length;
    auto [log_likelihood, log_likelihood_derivative] =
        this->LogLikelihoodAndDerivative(op);
    return std::make_pair(-log_likelihood, -branch_length * log_likelihood_derivative);
  };

  double current_log_branch_length = log(branch_lengths_(op.gpcsp_));
  double current_neg_log_likelihood =
      negative_log_likelihood_and_derivative(current_log_branch_length).first;
  const auto [log_branch_length, neg_log_likelihood] =
      Optimization::BrentMinimize<true>(
          negative_log_likelihood_and_derivative, current_log_branch_length,
          min_log_branch_length_, max_log_branch_length_,
          significant_digits_for_optimization_, max_iter_for_optimization_,
          step_size_for_log_space_optimization_);

  if (neg_log_likelihood > current_neg_log_likelihood) {
    branch_lengths_(op.gpcsp_) = exp(current_log_branch_length);
  } else {
    branch_lengths_(op.gpcsp_) = exp(log_branch_length);
  }
  branch_length_differences_(op.gpcsp_) =
      abs(exp(current_log_branch_length) - branch_lengths_(op.gpcsp_));
}

void GPEngine::GradientAscentOptimization(
    const GPOperations::OptimizeBranchLength& op) {
  auto log_likelihood_and_derivative = [this, &op](double branch_length) {
    branch_lengths_(op.gpcsp_) = branch_length;
    return this->LogLikelihoodAndDerivative(op);
  };
  const auto branch_length = Optimization::GradientAscent(
      log_likelihood_and_derivative, branch_lengths_(op.gpcsp_),
      significant_digits_for_optimization_, step_size_for_optimization_,
      min_log_branch_length_, max_iter_for_optimization_);
  branch_lengths_(op.gpcsp_) = branch_length;
}

void GPEngine::LogSpaceGradientAscentOptimization(
    const GPOperations::OptimizeBranchLength& op) {
  auto log_likelihood_and_derivative = [this, &op](double branch_length) {
    branch_lengths_(op.gpcsp_) = branch_length;
    return this->LogLikelihoodAndDerivative(op);
  };
  const auto branch_length = Optimization::LogSpaceGradientAscent(
      log_likelihood_and_derivative, branch_lengths_(op.gpcsp_),
      significant_digits_for_optimization_, step_size_for_log_space_optimization_,
      exp(min_log_branch_length_), max_iter_for_optimization_);
  branch_lengths_(op.gpcsp_) = branch_length;
}

void GPEngine::NewtonOptimization(const GPOperations::OptimizeBranchLength& op) {
  if (branch_length_differences_(op.gpcsp_) < branch_length_difference_threshold_) {
    return;
  }

  auto log_likelihood_and_first_two_derivatives = [this,
                                                   &op](double log_branch_length) {
    double x = exp(log_branch_length);
    branch_lengths_(op.gpcsp_) = x;
    auto [f_x, f_prime_x, f_double_prime_x] =
        this->LogLikelihoodAndFirstTwoDerivatives(op);
    // x = exp(y) --> f'(exp(y)) = exp(y) * f'(exp(y)) = x * f'(x)
    double f_prime_y = x * f_prime_x;
    double f_double_prime_y = f_prime_y + std::pow(x, 2) * f_double_prime_x;
    return std::make_tuple(f_x, f_prime_y, f_double_prime_y);
  };

  double current_log_branch_length = log(branch_lengths_(op.gpcsp_));
  const auto log_branch_length = Optimization::NewtonRaphsonOptimization(
      log_likelihood_and_first_two_derivatives, current_log_branch_length,
      significant_digits_for_optimization_, denominator_tolerance_for_newton_,
      min_log_branch_length_, max_log_branch_length_, max_iter_for_optimization_);

  branch_lengths_(op.gpcsp_) = exp(log_branch_length);

  branch_length_differences_(op.gpcsp_) =
      abs(exp(current_log_branch_length) - branch_lengths_(op.gpcsp_));
}

void GPEngine::HotStartBranchLengths(const RootedTreeCollection& tree_collection,
                                     const BitsetSizeMap& indexer) {
  size_t unique_gpcsp_count = branch_lengths_.size();
  branch_lengths_.setZero();

  EigenVectorXi observed_gpcsp_counts = EigenVectorXi::Zero(unique_gpcsp_count);
  // Set the branch length vector to be the total of the branch lengths for each PCSP,
  // and count the number of times we have seen each PCSP (into gpcsp_counts).
  auto tally_branch_lengths_and_gpcsp_count =
      [&observed_gpcsp_counts, this](size_t gpcsp_idx, const RootedTree& tree,
                                     const Node* focal_node) {
        branch_lengths_(gpcsp_idx) += tree.BranchLength(focal_node);
        observed_gpcsp_counts(gpcsp_idx)++;
      };
  FunctionOverRootedTreeCollection(tally_branch_lengths_and_gpcsp_count,
                                   tree_collection, indexer);
  for (size_t gpcsp_idx = 0; gpcsp_idx < unique_gpcsp_count; gpcsp_idx++) {
    if (observed_gpcsp_counts(gpcsp_idx) == 0) {
      branch_lengths_(gpcsp_idx) = default_branch_length_;
    } else {
      // Normalize the branch length total using the counts to get a mean branch
      // length.
      branch_lengths_(gpcsp_idx) /=
          static_cast<double>(observed_gpcsp_counts(gpcsp_idx));
    }
  }
}

SizeDoubleVectorMap GPEngine::GatherBranchLengths(
    const RootedTreeCollection& tree_collection, const BitsetSizeMap& indexer) {
  SizeDoubleVectorMap gpcsp_branchlengths_map;
  auto gather_branch_lengths = [&gpcsp_branchlengths_map](size_t gpcsp_idx,
                                                          const RootedTree& tree,
                                                          const Node* focal_node) {
    gpcsp_branchlengths_map[gpcsp_idx].push_back(tree.BranchLength(focal_node));
  };
  FunctionOverRootedTreeCollection(gather_branch_lengths, tree_collection, indexer);
  return gpcsp_branchlengths_map;
}

void GPEngine::FunctionOverRootedTreeCollection(
    std::function<void(size_t, const RootedTree&, const Node*)>
        function_on_tree_node_by_gpcsp,
    const RootedTreeCollection& tree_collection, const BitsetSizeMap& indexer) {
  const auto leaf_count = tree_collection.TaxonCount();
  const size_t default_index = branch_lengths_.size();
  for (const auto& tree : tree_collection.Trees()) {
    tree.Topology()->RootedPCSPPreorder(
        [&leaf_count, &default_index, &indexer, &tree, &function_on_tree_node_by_gpcsp](
            const Node* sister_node, const Node* focal_node, const Node* child0_node,
            const Node* child1_node) {
          Bitset gpcsp_bitset =
              SBNMaps::PCSPBitsetOf(leaf_count, sister_node, false, focal_node, false,
                                    child0_node, false, child1_node, false);
          const auto gpcsp_idx = AtWithDefault(indexer, gpcsp_bitset, default_index);
          if (gpcsp_idx != default_index) {
            function_on_tree_node_by_gpcsp(gpcsp_idx, tree, focal_node);
          }
        },
        true);
  }
}

void GPEngine::TakeFirstBranchLength(const RootedTreeCollection& tree_collection,
                                     const BitsetSizeMap& indexer) {
  size_t unique_gpcsp_count = branch_lengths_.size();
  branch_lengths_.setZero();
  EigenVectorXi observed_gpcsp_counts = EigenVectorXi::Zero(unique_gpcsp_count);
  // Set the branch length vector to be the first encountered branch length for each
  // PCSP, and mark when we have seen each PCSP (into observed_gpcsp_counts).
  auto set_first_branch_length_and_increment_gpcsp_count =
      [&observed_gpcsp_counts, this](size_t gpcsp_idx, const RootedTree& tree,
                                     const Node* focal_node) {
        if (observed_gpcsp_counts(gpcsp_idx) == 0) {
          branch_lengths_(gpcsp_idx) = tree.BranchLength(focal_node);
          observed_gpcsp_counts(gpcsp_idx)++;
        }
      };
  FunctionOverRootedTreeCollection(set_first_branch_length_and_increment_gpcsp_count,
                                   tree_collection, indexer);
  // If a branch length was not set above, set it to the default length.
  for (size_t gpcsp_idx = 0; gpcsp_idx < unique_gpcsp_count; gpcsp_idx++) {
    if (observed_gpcsp_counts(gpcsp_idx) == 0) {
      branch_lengths_(gpcsp_idx) = default_branch_length_;
    }
  }
}

EigenVectorXd GPEngine::CalculateQuartetHybridLikelihoods(
    const QuartetHybridRequest& request) {
  auto CheckRescaling = [this](size_t plv_idx) {
    Assert(rescaling_counts_[plv_idx] == 0,
           "Rescaling not implemented in CalculateQuartetHybridLikelihoods.");
  };
  std::vector<double> result;
  for (const auto& rootward_tip : request.rootward_tips_) {
    CheckRescaling(rootward_tip.plv_idx_);
    const double rootward_tip_prior =
        unconditional_node_probabilities_[rootward_tip.tip_node_id_];
    const double log_rootward_tip_prior = log(rootward_tip_prior);
    // #328 note that for the general case we should transpose the transition matrix
    // when coming down the tree.
    SetTransitionMatrixToHaveBranchLength(branch_lengths_[rootward_tip.gpcsp_idx_]);
    quartet_root_plv_ = transition_matrix_ * plv_handler_.GetPLV(rootward_tip.plv_idx_);
    for (const auto& sister_tip : request.sister_tips_) {
      CheckRescaling(sister_tip.plv_idx_);
      // Form the PLV on the root side of the central edge.
      SetTransitionMatrixToHaveBranchLength(branch_lengths_[sister_tip.gpcsp_idx_]);
      quartet_r_s_plv_.array() =
          quartet_root_plv_.array() *
          (transition_matrix_ * plv_handler_.GetPLV(sister_tip.plv_idx_)).array();
      // Advance it along the edge.
      SetTransitionMatrixToHaveBranchLength(
          branch_lengths_[request.central_gpcsp_idx_]);
      quartet_q_s_plv_ = transition_matrix_ * quartet_r_s_plv_;
      for (const auto& rotated_tip : request.rotated_tips_) {
        CheckRescaling(rotated_tip.plv_idx_);
        // Form the PLV on the root side of the sorted edge.
        SetTransitionMatrixToHaveBranchLength(branch_lengths_[rotated_tip.gpcsp_idx_]);
        quartet_r_sorted_plv_.array() =
            quartet_q_s_plv_.array() *
            (transition_matrix_ * plv_handler_.GetPLV(rotated_tip.plv_idx_)).array();
        for (const auto& sorted_tip : request.sorted_tips_) {
          CheckRescaling(sorted_tip.plv_idx_);
          // P(sigma_{ijkl} | \eta)
          const double non_sequence_based_log_probability = log(
              inverted_sbn_prior_[rootward_tip.gpcsp_idx_] * q_[sister_tip.gpcsp_idx_] *
              q_[rotated_tip.gpcsp_idx_] * q_[sorted_tip.gpcsp_idx_]);
          // Now calculate the sequence-based likelihood.
          SetTransitionMatrixToHaveBranchLength(branch_lengths_[sorted_tip.gpcsp_idx_]);
          per_pattern_log_likelihoods_ =
              (quartet_r_sorted_plv_.transpose() * transition_matrix_ *
               plv_handler_.GetPLV(sorted_tip.plv_idx_))
                  .diagonal()
                  .array()
                  .log();
          per_pattern_log_likelihoods_.array() -= log_rootward_tip_prior;
          result.push_back(non_sequence_based_log_probability +
                           per_pattern_log_likelihoods_.dot(site_pattern_weights_));
        }
      }
    }
  }
  return EigenVectorXdOfStdVectorDouble(result);
}

void GPEngine::ProcessQuartetHybridRequest(const QuartetHybridRequest& request) {
  if (request.IsFullyFormed()) {
    EigenVectorXd hybrid_log_likelihoods = CalculateQuartetHybridLikelihoods(request);
    hybrid_marginal_log_likelihoods_[request.central_gpcsp_idx_] =
        NumericalUtils::LogSum(hybrid_log_likelihoods);
  }
}

// ** I/O

std::string GPEngine::PLVToString(size_t plv_idx) const {
  std::stringstream out;
  for (auto&& row : plv_handler_.GetPLV(plv_idx).rowwise()) {
    out << row << std::endl;
  }
  return out.str();
}

void GPEngine::PrintPLV(size_t plv_idx) const {
  std::cout << PLVToString(plv_idx) << std::endl;
}
