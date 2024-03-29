// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"

#include "sugar.hpp"
#include "sbn_maps.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t node_count, size_t gpcsp_count,
                   const std::string& mmap_file_path, double rescaling_threshold,
                   EigenVectorXd sbn_prior,
                   EigenVectorXd unconditional_node_probabilities,
                   EigenVectorXd inverted_sbn_prior, bool use_gradients)
    : site_pattern_(std::move(site_pattern)),
      rescaling_threshold_(rescaling_threshold),
      log_rescaling_threshold_(log(rescaling_threshold)),
      plv_handler_(mmap_file_path, 0, site_pattern_.PatternCount(), resizing_factor_),
      unconditional_node_probabilities_(std::move(unconditional_node_probabilities)),
      branch_handler_(0),
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
  branch_handler_.SetCount(GetGPCSPCount());
  branch_handler_.SetSpareCount(GetSpareGPCSPCount());
  GrowGPCSPs(gpcsp_count, std::nullopt, std::nullopt, true);
  // Initialize PLV temporaries.
  quartet_root_plv_ = GetPLV(PVId(0));
  quartet_root_plv_.setZero();
  quartet_r_s_plv_ = quartet_root_plv_;
  quartet_q_s_plv_ = quartet_root_plv_;
  quartet_r_sorted_plv_ = quartet_root_plv_;

  InitializeBranchLengthHandler();
  UseGradientOptimization(use_gradients);
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
                        std::optional<const size_t> explicit_alloc,
                        const bool on_init) {
  const size_t old_node_count = GetNodeCount();
  const size_t old_plv_count = GetPLVCount();
  SetNodeCount(new_node_count);
  // Reallocate more space if needed.
  if ((GetPaddedNodeCount() > GetAllocatedNodeCount()) || explicit_alloc.has_value()) {
    SetAllocatedNodeCount(
        size_t(ceil(double(GetPaddedNodeCount()) * resizing_factor_)));
    if (explicit_alloc.has_value()) {
      Assert(explicit_alloc.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      SetAllocatedNodeCount(explicit_alloc.value() + GetSpareNodeCount());
    }
    plv_handler_.Resize(new_node_count, GetAllocatedNodeCount());
    rescaling_counts_.conservativeResize(GetAllocatedPLVCount());
    unconditional_node_probabilities_.conservativeResize(GetAllocatedNodeCount());
  }
  // Resize to fit without deallocating unused memory.
  rescaling_counts_.conservativeResize(GetPaddedPLVCount());
  if (on_init) {
    unconditional_node_probabilities_.conservativeResize(GetNodeCount());
  } else {
    unconditional_node_probabilities_.conservativeResize(GetPaddedNodeCount());
  }
  // Initialize new work space.
  Assert((GetPLVs().back().rows() == MmappedNucleotidePLV::base_count_) &&
             (GetPLVs().back().cols() ==
              static_cast<Eigen::Index>(site_pattern_.PatternCount())) &&
             (size_t(GetPLVs().size()) == GetAllocatedPLVCount()),
         "Didn't get the right shape of PLVs out of Subdivide.");
  for (PVId pv_id = PVId(old_plv_count); pv_id < GetPaddedPLVCount(); pv_id++) {
    rescaling_counts_[pv_id.value_] = 0;
    GetPLV(pv_id).setZero();
  }
  for (NodeId node_id = NodeId(old_node_count); node_id < GetNodeCount(); node_id++) {
    if (!on_init) {
      unconditional_node_probabilities_[node_id.value_] = 1.;
    }
  }
  // Reindex work space to realign with DAG.
  if (node_reindexer.has_value()) {
    ReindexPLVs(node_reindexer.value(), old_node_count);
  }
}

void GPEngine::GrowGPCSPs(const size_t new_gpcsp_count,
                          std::optional<const Reindexer> gpcsp_reindexer,
                          std::optional<const size_t> explicit_alloc,
                          const bool on_init) {
  const size_t old_gpcsp_count = GetGPCSPCount();
  SetGPCSPCount(new_gpcsp_count);
  branch_handler_.Resize(new_gpcsp_count, std::nullopt, explicit_alloc,
                         gpcsp_reindexer);
  // Reallocate more space if needed.
  if ((GetPaddedGPCSPCount() > GetAllocatedGPCSPCount()) ||
      explicit_alloc.has_value()) {
    SetAllocatedGPCSPCount(
        size_t(ceil(double(GetPaddedGPCSPCount()) * resizing_factor_)));
    if (explicit_alloc.has_value()) {
      Assert(explicit_alloc.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      SetAllocatedGPCSPCount(explicit_alloc.value() + GetSpareGPCSPCount());
    }
    hybrid_marginal_log_likelihoods_.conservativeResize(GetAllocatedGPCSPCount());
    log_likelihoods_.conservativeResize(GetAllocatedGPCSPCount(),
                                        site_pattern_.PatternCount());
    q_.conservativeResize(GetAllocatedGPCSPCount());
    inverted_sbn_prior_.conservativeResize(GetAllocatedGPCSPCount());
  }
  // Resize to fit without deallocating unused memory.
  hybrid_marginal_log_likelihoods_.conservativeResize(GetPaddedGPCSPCount());
  log_likelihoods_.conservativeResize(GetPaddedGPCSPCount(),
                                      site_pattern_.PatternCount());
  if (on_init) {
    q_.conservativeResize(GetGPCSPCount());
    inverted_sbn_prior_.conservativeResize(GetGPCSPCount());
  } else {
    q_.conservativeResize(GetPaddedGPCSPCount());
    inverted_sbn_prior_.conservativeResize(GetPaddedGPCSPCount());
  }
  // Initialize new work space.
  for (size_t i = old_gpcsp_count; i < GetPaddedGPCSPCount(); i++) {
    hybrid_marginal_log_likelihoods_[i] = DOUBLE_NEG_INF;
  }
  if (on_init) {
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
      plv_handler_.BuildPVReindexer(node_reindexer, old_node_count, GetNodeCount());
  // Reindex data vectors
  plv_handler_.Reindex(plv_reindexer);
  Reindexer::ReindexInPlace<EigenVectorXi, int>(rescaling_counts_, plv_reindexer,
                                                GetPLVCount());
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
  Reindexer::ReindexInPlace<EigenVectorXd, double>(hybrid_marginal_log_likelihoods_,
                                                   gpcsp_reindexer, GetGPCSPCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(q_, gpcsp_reindexer,
                                                   GetGPCSPCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(inverted_sbn_prior_, gpcsp_reindexer,
                                                   GetGPCSPCount());
}

void GPEngine::GrowSparePLVs(const size_t new_node_spare_count) {
  if (new_node_spare_count > GetSpareNodeCount()) {
    SetSpareNodeCount(new_node_spare_count);
    GrowPLVs(GetNodeCount());
  }
}

void GPEngine::GrowSpareGPCSPs(const size_t new_gpcsp_spare_count) {
  if (new_gpcsp_spare_count > GetSpareGPCSPCount()) {
    SetSpareGPCSPCount(new_gpcsp_spare_count);
    branch_handler_.SetSpareCount(new_gpcsp_spare_count);
    GrowGPCSPs(GetGPCSPCount());
  }
}

// ** GPOperations

void GPEngine::operator()(const GPOperations::ZeroPLV& op) {
  GetPLV(PVId(op.dest_)).setZero();
  rescaling_counts_(op.dest_) = 0;
}

void GPEngine::operator()(const GPOperations::SetToStationaryDistribution& op) {
  auto& plv = GetPLV(PVId(op.dest_));
  for (Eigen::Index row_idx = 0; row_idx < plv.rows(); ++row_idx) {
    // Multiplication by q_ avoids special treatment of the rhat vector for the
    // rootsplits.
    plv.row(row_idx).array() =
        q_(op.root_gpcsp_idx_) * stationary_distribution_(row_idx);
  }
  rescaling_counts_(op.dest_) = 0;
}

void GPEngine::operator()(const GPOperations::IncrementWithWeightedEvolvedPLV& op) {
  const auto branch_length = branch_handler_(EdgeId(op.gpcsp_));
  SetTransitionMatrixToHaveBranchLength(branch_length);
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
  GetPLV(PVId(op.dest_)) +=
      rescaling_factor * q_(op.gpcsp_) * transition_matrix_ * GetPLV(PVId(op.src_));
}

void GPEngine::operator()(const GPOperations::ResetMarginalLikelihood& op) {  // NOLINT
  ResetLogMarginalLikelihood();
}

void GPEngine::operator()(const GPOperations::IncrementMarginalLikelihood& op) {
  Assert(rescaling_counts_(op.stationary_times_prior_) == 0,
         "Surprise! Rescaled stationary distribution in IncrementMarginalLikelihood");
  // This operation does two things: increment the overall per-site log marginal
  // likelihood, and also set the conditional per-rootsplit marginal likelihood.
  //
  // We first calculate the unconditional contribution of the rootsplit to the overall
  // per-site marginal likelihood. It's an unconditional contribution because our
  // stationary distribution incorporates the prior on rootsplits.
  log_likelihoods_.row(op.rootsplit_) =
      (GetPLV(PVId(op.stationary_times_prior_)).transpose() * GetPLV(PVId(op.p_)))
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
  GetPLV(PVId(op.dest_)).array() =
      GetPLV(PVId(op.src1_)).array() * GetPLV(PVId(op.src2_)).array();
  rescaling_counts_(op.dest_) =
      rescaling_counts_(op.src1_) + rescaling_counts_(op.src2_);
  AssertPLVIsFinite(op.dest_, "Multiply dest_ is not finite");
  RescalePLVIfNeeded(op.dest_);
}

void GPEngine::operator()(const GPOperations::Likelihood& op) {
  SetTransitionMatrixToHaveBranchLength(branch_handler_(EdgeId(op.dest_)));
  PreparePerPatternLogLikelihoodsForGPCSP(op.parent_, op.child_);
  log_likelihoods_.row(op.dest_) = per_pattern_log_likelihoods_;
}

void GPEngine::operator()(const GPOperations::OptimizeBranchLength& op) {
  return OptimizeBranchLength(op);
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
  branch_handler_.GetBranchLengths().GetData().segment(0, GetGPCSPCount()) =
      branch_lengths;
};

void GPEngine::SetBranchLengthsToConstant(double branch_length) {
  branch_handler_.GetBranchLengths().GetData().setConstant(branch_length);
};

void GPEngine::SetBranchLengthsToDefault() {
  branch_handler_.GetBranchLengths().GetData().setConstant(
      branch_handler_.GetDefaultBranchLength());
};

void GPEngine::ResetLogMarginalLikelihood() {
  log_marginal_likelihood_.setConstant(DOUBLE_NEG_INF);
}

void GPEngine::CopyNodeData(const NodeId src_node_idx, const NodeId dest_node_idx) {
  Assert(
      (src_node_idx < GetPaddedNodeCount()) && (dest_node_idx < GetPaddedNodeCount()),
      "Cannot copy node data with src or dest index out-of-range.");
  unconditional_node_probabilities_[dest_node_idx.value_] =
      unconditional_node_probabilities_[src_node_idx.value_];
}

void GPEngine::CopyPLVData(const size_t src_plv_idx, const size_t dest_plv_idx) {
  Assert((src_plv_idx < GetPaddedPLVCount()) && (dest_plv_idx < GetPaddedPLVCount()),
         "Cannot copy PLV data with src or dest index out-of-range.");
  GetPLV(PVId(dest_plv_idx)) = GetPLV(PVId(src_plv_idx));
  rescaling_counts_[dest_plv_idx] = rescaling_counts_[src_plv_idx];
}

void GPEngine::CopyGPCSPData(const EdgeId src_gpcsp_idx, const EdgeId dest_gpcsp_idx) {
  Assert((src_gpcsp_idx < GetPaddedGPCSPCount()) &&
             (dest_gpcsp_idx < GetPaddedGPCSPCount()),
         "Cannot copy PLV data with src or dest index out-of-range.");
  branch_handler_(dest_gpcsp_idx) = branch_handler_(src_gpcsp_idx);
  q_[dest_gpcsp_idx.value_] = q_[src_gpcsp_idx.value_];
  inverted_sbn_prior_[dest_gpcsp_idx.value_] =
      inverted_sbn_prior_[src_gpcsp_idx.value_];
}

// ** Access

double GPEngine::GetLogMarginalLikelihood() const {
  return (log_marginal_likelihood_.array() * site_pattern_weights_.array()).sum();
}

EigenVectorXd GPEngine::GetBranchLengths() const {
  return branch_handler_.GetBranchLengths().GetData().segment(0, GetGPCSPCount());
};

EigenVectorXd GPEngine::GetBranchLengths(const size_t start,
                                         const size_t length) const {
  Assert(start + length <= GetPaddedGPCSPCount(),
         "Requested range of BranchLengths is out-of-range.");
  return branch_handler_.GetBranchLengths().GetData().segment(start, length);
};

EigenVectorXd GPEngine::GetSpareBranchLengths(const size_t start,
                                              const size_t length) const {
  return GetBranchLengths(GetSpareGPCSPIndex(start), length);
}

EigenVectorXd GPEngine::GetBranchLengthDifferences() const {
  return branch_handler_.GetBranchDifferences().GetData().segment(0, GetGPCSPCount());
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

EigenVectorXd GPEngine::GetSparePerGPCSPLogLikelihoods(const size_t start,
                                                       const size_t length) const {
  return GetPerGPCSPLogLikelihoods(GetSpareGPCSPIndex(start), length);
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

DoublePair GPEngine::LogLikelihoodAndDerivative(
    const GPOperations::OptimizeBranchLength& op) {
  return LogLikelihoodAndDerivative(op.gpcsp_, op.rootward_, op.leafward_);
}

DoublePair GPEngine::LogLikelihoodAndDerivative(const size_t gpcsp,
                                                const size_t rootward,
                                                const size_t leafward) {
  SetTransitionAndDerivativeMatricesToHaveBranchLength(branch_handler_(EdgeId(gpcsp)));
  PreparePerPatternLogLikelihoodsForGPCSP(rootward, leafward);
  // The prior is expressed using the current value of q_.
  // The phylogenetic component of the likelihood is weighted with the number of times
  // we see the site patterns.
  const double log_likelihood = per_pattern_log_likelihoods_.dot(site_pattern_weights_);

  // The per-site likelihood derivative is calculated in the same way as the per-site
  // likelihood, but using the derivative matrix instead of the transition matrix.
  // We first prepare two useful vectors _without_ likelihood rescaling, because the
  // rescalings cancel out in the ratio below.
  PrepareUnrescaledPerPatternLikelihoodDerivatives(rootward, leafward);
  PrepareUnrescaledPerPatternLikelihoods(rootward, leafward);
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
  return LogLikelihoodAndFirstTwoDerivatives(op.gpcsp_, op.rootward_, op.leafward_);
}

std::tuple<double, double, double> GPEngine::LogLikelihoodAndFirstTwoDerivatives(
    const size_t gpcsp, const size_t rootward, const size_t leafward) {
  SetTransitionAndDerivativeMatricesToHaveBranchLength(branch_handler_(EdgeId(gpcsp)));
  PreparePerPatternLogLikelihoodsForGPCSP(rootward, leafward);

  const double log_likelihood = per_pattern_log_likelihoods_.dot(site_pattern_weights_);

  // The per-site likelihood derivative is calculated in the same way as the per-site
  // likelihood, but using the derivative matrix instead of the transition matrix.
  // We first prepare two useful vectors _without_ likelihood rescaling, because the
  // rescalings cancel out in the ratio below.
  PrepareUnrescaledPerPatternLikelihoodDerivatives(rootward, leafward);
  PrepareUnrescaledPerPatternLikelihoods(rootward, leafward);
  // If l_i is the per-site likelihood, the derivative of log(l_i) is the derivative
  // of l_i divided by l_i.
  per_pattern_likelihood_derivative_ratios_ =
      per_pattern_likelihood_derivatives_.array() / per_pattern_likelihoods_.array();
  const double log_likelihood_gradient =
      per_pattern_likelihood_derivative_ratios_.dot(site_pattern_weights_);

  // Second derivative is calculated the same way, but has an extra term due to
  // the product rule.

  PrepareUnrescaledPerPatternLikelihoodSecondDerivatives(rootward, leafward);

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
  for (auto& plv : GetPLVs()) {
    plv.setZero();
  }
  NodeId taxon_idx = 0;
  for (const auto& pattern : site_pattern_.GetPatterns()) {
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative symbol!");
      if (symbol == MmappedNucleotidePLV::base_count_) {  // Gap character.
        GetPLV(PVId(taxon_idx.value_)).col(site_idx).setConstant(1.);
      } else if (symbol < MmappedNucleotidePLV::base_count_) {
        GetPLV(PVId(taxon_idx.value_))(symbol, site_idx) = 1.;
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
  GetPLV(PVId(plv_idx)) /=
      pow(rescaling_threshold_, static_cast<double>(rescaling_count));
  rescaling_counts_(plv_idx) += rescaling_count;
}

void GPEngine::AssertPLVIsFinite(size_t plv_idx, const std::string& message) const {
  Assert(GetPLV(PVId(plv_idx)).array().isFinite().all(), message);
}

std::pair<double, double> GPEngine::PLVMinMax(size_t plv_idx) const {
  return {GetPLV(PVId(plv_idx)).minCoeff(), GetPLV(PVId(plv_idx)).maxCoeff()};
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

void GPEngine::InitializeBranchLengthHandler() {
  // Set Nongradient Brent.
  DAGBranchHandler::NegLogLikelihoodFunc brent_nongrad_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double log_branch_length) {
        SetTransitionMatrixToHaveBranchLength(exp(log_branch_length));
        PreparePerPatternLogLikelihoodsForGPCSP(parent_id.value_, child_id.value_);
        return -per_pattern_log_likelihoods_.dot(site_pattern_weights_);
      };
  branch_handler_.SetBrentFunc(brent_nongrad_func);
  // Set Gradient Brent.
  DAGBranchHandler::NegLogLikelihoodAndDerivativeFunc brent_grad_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double log_branch_length) {
        double branch_length = exp(log_branch_length);
        branch_handler_(edge_id) = branch_length;
        auto [log_likelihood, log_likelihood_derivative] =
            this->LogLikelihoodAndDerivative(edge_id.value_, parent_id.value_,
                                             child_id.value_);
        return std::make_pair(-log_likelihood,
                              -branch_length * log_likelihood_derivative);
      };
  branch_handler_.SetBrentWithGradientFunc(brent_grad_func);
  // Set Gradient Ascent.
  DAGBranchHandler::LogLikelihoodAndDerivativeFunc grad_ascent_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double branch_length) {
        branch_handler_(edge_id) = branch_length;
        return this->LogLikelihoodAndDerivative(edge_id.value_, parent_id.value_,
                                                child_id.value_);
      };
  branch_handler_.SetGradientAscentFunc(grad_ascent_func);
  // Set Logspace Gradient Ascent.
  DAGBranchHandler::LogLikelihoodAndDerivativeFunc logspace_grad_ascent_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double branch_length) {
        branch_handler_(edge_id) = branch_length;
        return this->LogLikelihoodAndDerivative(edge_id.value_, parent_id.value_,
                                                child_id.value_);
      };
  branch_handler_.SetLogSpaceGradientAscentFunc(logspace_grad_ascent_func);
  // Set Newton-Raphson.
  DAGBranchHandler::LogLikelihoodAndFirstTwoDerivativesFunc newton_raphson_func =
      [this](EdgeId edge_id, PVId parent_id, PVId child_id, double log_branch_length) {
        double x = exp(log_branch_length);
        branch_handler_(edge_id) = x;
        auto [f_x, f_prime_x, f_double_prime_x] =
            this->LogLikelihoodAndFirstTwoDerivatives(edge_id.value_, parent_id.value_,
                                                      child_id.value_);
        // x = exp(y) --> f'(exp(y)) = exp(y) * f'(exp(y)) = x * f'(x)
        double f_prime_y = x * f_prime_x;
        double f_double_prime_y = f_prime_y + std::pow(x, 2) * f_double_prime_x;
        return std::make_tuple(f_x, f_prime_y, f_double_prime_y);
      };
  branch_handler_.SetNewtonRaphsonFunc(newton_raphson_func);
}

void GPEngine::SetOptimizationMethod(const OptimizationMethod method) {
  branch_handler_.SetOptimizationMethod(method);
}

void GPEngine::UseGradientOptimization(const bool use_gradients) {
  auto optimization_method =
      (use_gradients ? OptimizationMethod::BrentOptimizationWithGradients
                     : OptimizationMethod::BrentOptimization);
  branch_handler_.SetOptimizationMethod(optimization_method);
}

void GPEngine::OptimizeBranchLength(const GPOperations::OptimizeBranchLength& op) {
  branch_handler_.OptimizeBranchLength(EdgeId(op.gpcsp_), PVId(op.rootward_),
                                       PVId(op.leafward_), !IsFirstOptimization());
}

void GPEngine::SetSignificantDigitsForOptimization(int significant_digits) {
  branch_handler_.SetSignificantDigitsForOptimization(significant_digits);
}

void GPEngine::HotStartBranchLengths(const RootedTreeCollection& tree_collection,
                                     const BitsetSizeMap& indexer) {
  size_t unique_gpcsp_count = branch_handler_.size();
  branch_handler_.GetBranchLengths().GetData().setZero();

  EigenVectorXi observed_gpcsp_counts = EigenVectorXi::Zero(unique_gpcsp_count);
  // Set the branch length vector to be the total of the branch lengths for each PCSP,
  // and count the number of times we have seen each PCSP (into gpcsp_counts).
  auto tally_branch_handler_and_gpcsp_count =
      [&observed_gpcsp_counts, this](EdgeId gpcsp_idx, const Bitset& bitset,
                                     const RootedTree& tree, const size_t tree_id,
                                     const Node* focal_node) {
        branch_handler_(gpcsp_idx) += tree.BranchLength(focal_node);
        observed_gpcsp_counts(gpcsp_idx.value_)++;
      };
  RootedSBNMaps::FunctionOverRootedTreeCollection(tally_branch_handler_and_gpcsp_count,
                                                  tree_collection, indexer,
                                                  branch_handler_.size());
  for (EdgeId gpcsp_idx = 0; gpcsp_idx.value_ < unique_gpcsp_count;
       gpcsp_idx.value_++) {
    if (observed_gpcsp_counts(gpcsp_idx.value_) == 0) {
      branch_handler_(gpcsp_idx) = branch_handler_.GetDefaultBranchLength();
    } else {
      // Normalize the branch length total using the counts to get a mean branch
      // length.
      branch_handler_(gpcsp_idx) /=
          static_cast<double>(observed_gpcsp_counts(gpcsp_idx.value_));
    }
  }
}

SizeDoubleVectorMap GPEngine::GatherBranchLengths(
    const RootedTreeCollection& tree_collection, const BitsetSizeMap& indexer) {
  SizeDoubleVectorMap gpcsp_branchlengths_map;
  auto gather_branch_lengths = [&gpcsp_branchlengths_map](
                                   EdgeId gpcsp_idx, const Bitset& bitset,
                                   const RootedTree& tree, const size_t tree_id,
                                   const Node* focal_node) {
    gpcsp_branchlengths_map[gpcsp_idx.value_].push_back(tree.BranchLength(focal_node));
  };
  RootedSBNMaps::FunctionOverRootedTreeCollection(
      gather_branch_lengths, tree_collection, indexer, branch_handler_.size());
  return gpcsp_branchlengths_map;
}

void GPEngine::TakeFirstBranchLength(const RootedTreeCollection& tree_collection,
                                     const BitsetSizeMap& indexer) {
  size_t unique_gpcsp_count = branch_handler_.size();
  branch_handler_.GetBranchLengths().GetData().setZero();
  EigenVectorXi observed_gpcsp_counts = EigenVectorXi::Zero(unique_gpcsp_count);
  // Set the branch length vector to be the first encountered branch length for each
  // PCSP, and mark when we have seen each PCSP (into observed_gpcsp_counts).
  auto set_first_branch_length_and_increment_gpcsp_count =
      [&observed_gpcsp_counts, this](EdgeId gpcsp_idx, const Bitset& bitset,
                                     const RootedTree& tree, const size_t tree_id,
                                     const Node* focal_node) {
        if (observed_gpcsp_counts(gpcsp_idx.value_) == 0) {
          branch_handler_(gpcsp_idx) = tree.BranchLength(focal_node);
          observed_gpcsp_counts(gpcsp_idx.value_)++;
        }
      };
  RootedSBNMaps::FunctionOverRootedTreeCollection(
      set_first_branch_length_and_increment_gpcsp_count, tree_collection, indexer,
      branch_handler_.size());
  // If a branch length was not set above, set it to the default length.
  for (EdgeId gpcsp_idx = 0; gpcsp_idx < unique_gpcsp_count; gpcsp_idx.value_++) {
    if (observed_gpcsp_counts(gpcsp_idx.value_) == 0) {
      branch_handler_(gpcsp_idx) = branch_handler_.GetDefaultBranchLength();
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
    SetTransitionMatrixToHaveBranchLength(
        branch_handler_(EdgeId(rootward_tip.gpcsp_idx_)));
    quartet_root_plv_ = transition_matrix_ * GetPLV(PVId(rootward_tip.plv_idx_));
    for (const auto& sister_tip : request.sister_tips_) {
      CheckRescaling(sister_tip.plv_idx_);
      // Form the PLV on the root side of the central edge.
      SetTransitionMatrixToHaveBranchLength(
          branch_handler_(EdgeId(sister_tip.gpcsp_idx_)));
      quartet_r_s_plv_.array() =
          quartet_root_plv_.array() *
          (transition_matrix_ * GetPLV(PVId(sister_tip.plv_idx_))).array();
      // Advance it along the edge.
      SetTransitionMatrixToHaveBranchLength(
          branch_handler_(EdgeId(request.central_gpcsp_idx_)));
      quartet_q_s_plv_ = transition_matrix_ * quartet_r_s_plv_;
      for (const auto& rotated_tip : request.rotated_tips_) {
        CheckRescaling(rotated_tip.plv_idx_);
        // Form the PLV on the root side of the sorted edge.
        SetTransitionMatrixToHaveBranchLength(
            branch_handler_(EdgeId(rotated_tip.gpcsp_idx_)));
        quartet_r_sorted_plv_.array() =
            quartet_q_s_plv_.array() *
            (transition_matrix_ * GetPLV(PVId(rotated_tip.plv_idx_))).array();
        for (const auto& sorted_tip : request.sorted_tips_) {
          CheckRescaling(sorted_tip.plv_idx_);
          // P(sigma_{ijkl} | \eta)
          const double non_sequence_based_log_probability = log(
              inverted_sbn_prior_[rootward_tip.gpcsp_idx_] * q_[sister_tip.gpcsp_idx_] *
              q_[rotated_tip.gpcsp_idx_] * q_[sorted_tip.gpcsp_idx_]);
          // Now calculate the sequence-based likelihood.
          SetTransitionMatrixToHaveBranchLength(
              branch_handler_(EdgeId(sorted_tip.gpcsp_idx_)));
          per_pattern_log_likelihoods_ =
              (quartet_r_sorted_plv_.transpose() * transition_matrix_ *
               GetPLV(PVId(sorted_tip.plv_idx_)))
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

std::string GPEngine::PLVToString(const PVId plv_idx) const {
  return plv_handler_.ToString(plv_idx);
}

std::string GPEngine::LogLikelihoodMatrixToString() const {
  std::stringstream out;
  for (Eigen::Index i = 0; i < log_likelihoods_.rows(); i++) {
    for (Eigen::Index j = 0; j < log_likelihoods_.cols(); j++) {
      out << "[" << i << "," << j << "]: " << log_likelihoods_(i, j) << "\t";
    }
    out << std::endl;
  }
  return out.str();
}
