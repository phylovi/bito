// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"

#include "plv_handler.hpp"
#include "optimization.hpp"
#include "sugar.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t node_count,
                   size_t plv_count_per_node, size_t gpcsp_count,
                   const std::string& mmap_file_path, double rescaling_threshold,
                   EigenVectorXd sbn_prior,
                   EigenVectorXd unconditional_node_probabilities,
                   EigenVectorXd inverted_sbn_prior)
    : site_pattern_(std::move(site_pattern)),
      rescaling_threshold_(rescaling_threshold),
      log_rescaling_threshold_(log(rescaling_threshold)),
      plv_count_per_node_(plv_count_per_node),
      mmap_file_path_(mmap_file_path),
      mmapped_master_plv_(mmap_file_path_,
                          (node_count + node_padding_) * plv_count_per_node_ *
                              size_t(resizing_factor_) * site_pattern_.PatternCount()),
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
  quartet_root_plv_ = plvs_.at(0);
  quartet_root_plv_.setZero();
  quartet_r_s_plv_ = quartet_root_plv_;
  quartet_q_s_plv_ = quartet_root_plv_;
  quartet_r_sorted_plv_ = quartet_root_plv_;
}

void GPEngine::InitializePriors(EigenVectorXd sbn_prior,
                                EigenVectorXd unconditional_node_probabilities,
                                EigenVectorXd inverted_sbn_prior) {
  Assert(size_t(unconditional_node_probabilities.size()) == node_count_,
         "unconditional_node_probabilities is wrong size for GPEngine.");
  Assert(size_t(sbn_prior.size()) == gpcsp_count_,
         "sbn_prior is wrong size for GPEngine.");
  Assert(size_t(inverted_sbn_prior.size()) == gpcsp_count_,
         "inverted_sbn_prior is wrong size for GPEngine.");
  unconditional_node_probabilities_.segment(0, node_count_) =
      unconditional_node_probabilities;
  q_.segment(0, gpcsp_count_) = sbn_prior;
  inverted_sbn_prior_.segment(0, gpcsp_count_) = inverted_sbn_prior;
}

void GPEngine::GrowPLVs(const size_t new_node_count,
                        std::optional<const Reindexer> node_reindexer,
                        std::optional<const size_t> explicit_allocation,
                        const bool on_initialization) {
  const size_t old_node_count = GetNodeCount();
  const size_t old_plv_count = GetPLVCount();
  node_count_ = new_node_count;
  // Reallocate more space if needed.
  if ((GetPaddedNodeCount() > GetAllocatedNodeCount()) ||
      explicit_allocation.has_value()) {
    node_alloc_ = size_t(ceil(double(GetPaddedNodeCount()) * resizing_factor_));
    if (explicit_allocation.has_value()) {
      Assert(explicit_allocation.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      node_alloc_ = explicit_allocation.value() + node_padding_;
    }
    mmapped_master_plv_.Resize(GetAllocatedPLVCount() * site_pattern_.PatternCount());
    plvs_ = mmapped_master_plv_.Subdivide(GetAllocatedPLVCount());
    rescaling_counts_.conservativeResize(GetAllocatedPLVCount());
    unconditional_node_probabilities_.conservativeResize(GetAllocatedNodeCount());
  }
  // Resize to fit without deallocating unused memory.
  rescaling_counts_.conservativeResize(GetPaddedPLVCount());
  unconditional_node_probabilities_.conservativeResize(GetNodeCount());
  // Initialize new work space.
  Assert((plvs_.back().rows() == MmappedNucleotidePLV::base_count_) &&
             (plvs_.back().cols() ==
              static_cast<Eigen::Index>(site_pattern_.PatternCount())) &&
             (size_t(plvs_.size()) == GetAllocatedPLVCount()),
         "Didn't get the right shape of PLVs out of Subdivide.");
  for (size_t i = old_plv_count; i < GetPaddedPLVCount(); i++) {
    rescaling_counts_[i] = 0;
    plvs_.at(i).setZero();
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
  const size_t old_gpcsp_count = gpcsp_count_;
  gpcsp_count_ = new_gpcsp_count;
  // Reallocate more space if needed.
  if ((GetPaddedGPCSPCount() > gpcsp_alloc_) || explicit_allocation.has_value()) {
    gpcsp_alloc_ = size_t(ceil(double(GetPaddedGPCSPCount()) * resizing_factor_));
    if (explicit_allocation.has_value()) {
      Assert(explicit_allocation.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      gpcsp_alloc_ = explicit_allocation.value() + gpcsp_padding_;
    }
    branch_lengths_.conservativeResize(GetAllocatedGPCSPCount());
    hybrid_marginal_log_likelihoods_.conservativeResize(GetAllocatedGPCSPCount());
    log_likelihoods_.conservativeResize(GetAllocatedGPCSPCount(),
                                        site_pattern_.PatternCount());
    q_.conservativeResize(GetAllocatedGPCSPCount());
    inverted_sbn_prior_.conservativeResize(GetAllocatedGPCSPCount());
  }
  // Resize to fit without deallocating unused memory.
  branch_lengths_.conservativeResize(GetPaddedGPCSPCount());
  hybrid_marginal_log_likelihoods_.conservativeResize(GetPaddedGPCSPCount());
  log_likelihoods_.conservativeResize(GetPaddedGPCSPCount(),
                                      site_pattern_.PatternCount());
  q_.conservativeResize(GetGPCSPCount());
  inverted_sbn_prior_.conservativeResize(GetGPCSPCount());
  // Initialize new work space.
  for (size_t i = old_gpcsp_count; i < GetPaddedGPCSPCount(); i++) {
    branch_lengths_[i] = default_branch_length_;
    hybrid_marginal_log_likelihoods_[i] = DOUBLE_NEG_INF;
  }
  for (size_t i = old_gpcsp_count; i < GetGPCSPCount(); i++) {
    if (!on_initialization) {
      q_[i] = 1.;
      inverted_sbn_prior_[i] = 1.;
    }
  }
  // Reindex work space to realign with DAG.
  if (gpcsp_reindexer.has_value()) {
    ReindexGPCSPs(gpcsp_reindexer.value(), old_gpcsp_count);
  }
}

void GPEngine::ReindexPLVs(const Reindexer node_reindexer,
                           const size_t old_node_count) {
  Assert(node_reindexer.size() == node_count_,
         "Node Reindexer is the wrong size for GPEngine.");
  Assert(node_reindexer.IsValid(node_count_), "Node Reindexer is not valid.");

  // Expand node_reindexer into plv_reindexer.
  Reindexer plv_reindexer(node_count_ * plv_count_per_node_);
  size_t new_data_idx = old_node_count * plv_count_per_node_;
  for (size_t i = 0; i < node_count_; i++) {
    const size_t old_node_idx = i;
    const size_t new_node_idx = node_reindexer.GetNewIndexByOldIndex(old_node_idx);
    for (const auto plv_type : PLVHandler::PLVTypeIterator()) {
      // Either get input plv_index from old plvs, or get new plv_index (new data is
      // irrelevant, so just get next available index).
      size_t old_plv_idx;
      if (old_node_idx < old_node_count) {
        old_plv_idx = PLVHandler::GetPLVIndex(plv_type, old_node_count, old_node_idx);
      } else {
        old_plv_idx = new_data_idx;
        new_data_idx++;
      }
      const size_t new_plv_idx =
          PLVHandler::GetPLVIndex(plv_type, node_count_, new_node_idx);
      plv_reindexer.SetReindex(old_plv_idx, new_plv_idx);
    }
  }
  Assert(plv_reindexer.IsValid(GetPLVCount()), "PLV Reindexer is not valid.");
  // Reindex data vectors
  Reindexer::ReindexInPlace(plvs_, plv_reindexer, GetPLVCount(),
                            plvs_.at(GetPLVCount()), plvs_.at(GetPLVCount() + 1));
  Reindexer::ReindexInPlace<EigenVectorXi, int>(rescaling_counts_, plv_reindexer,
                                                GetPLVCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(unconditional_node_probabilities_,
                                                   node_reindexer, GetNodeCount());
}

void GPEngine::ReindexGPCSPs(const Reindexer gpcsp_reindexer,
                             const size_t old_gpcsp_count) {
  Assert(gpcsp_reindexer.size() == gpcsp_count_,
         "GPCSP Reindexer is the wrong size for GPEngine.");
  Assert(gpcsp_reindexer.IsValid(gpcsp_count_),
         "GPCSP Reindexer is not valid for GPEngine size.");
  // Reindex data vectors.
  Reindexer::ReindexInPlace<EigenVectorXd, double>(branch_lengths_, gpcsp_reindexer,
                                                   GetGPCSPCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(hybrid_marginal_log_likelihoods_,
                                                   gpcsp_reindexer, GetGPCSPCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(q_, gpcsp_reindexer,
                                                   GetGPCSPCount());
  Reindexer::ReindexInPlace<EigenVectorXd, double>(inverted_sbn_prior_, gpcsp_reindexer,
                                                   GetGPCSPCount());
}

void GPEngine::operator()(const GPOperations::ZeroPLV& op) {
  plvs_.at(op.dest_).setZero();
  rescaling_counts_(op.dest_) = 0;
}

void GPEngine::operator()(const GPOperations::SetToStationaryDistribution& op) {
  auto& plv = plvs_.at(op.dest_);
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
  // unavoidable without special-purpose truncation code, which doesn't seem worthwhile.
  plvs_.at(op.dest_) +=
      rescaling_factor * q_(op.gpcsp_) * transition_matrix_ * plvs_.at(op.src_);
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
      (plvs_.at(op.stationary_times_prior_).transpose() * plvs_.at(op.p_))
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
  plvs_.at(op.dest_).array() = plvs_.at(op.src1_).array() * plvs_.at(op.src2_).array();
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
  BrentOptimization(op);
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
  diagonal_matrix_.diagonal() = eigenvalues_.array() * diagonal_vector_.array();
  derivative_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void GPEngine::SetTransitionMatrixToHaveBranchLengthAndTranspose(double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ =
      inverse_eigenmatrix_.transpose() * diagonal_matrix_ * eigenmatrix_.transpose();
}

void GPEngine::SetBranchLengths(EigenVectorXd branch_lengths) {
  Assert(size_t(branch_lengths.size()) == gpcsp_count_,
         "Size mismatch in GPEngine::SetBranchLengths.");
  branch_lengths_.segment(0, gpcsp_count_) = branch_lengths;
};

void GPEngine::SetBranchLengthsToConstant(double branch_length) {
  branch_lengths_.setConstant(branch_length);
};

void GPEngine::ResetLogMarginalLikelihood() {
  log_marginal_likelihood_.setConstant(DOUBLE_NEG_INF);
}

double GPEngine::GetLogMarginalLikelihood() const {
  return (log_marginal_likelihood_.array() * site_pattern_weights_.array()).sum();
}

EigenVectorXd GPEngine::GetBranchLengths() const {
  return branch_lengths_.segment(0, gpcsp_count_);
};

EigenVectorXd GPEngine::GetPerGPCSPLogLikelihoods() const {
  return log_likelihoods_.block(0, 0, gpcsp_count_, log_likelihoods_.cols()) *
         site_pattern_weights_;
};

EigenVectorXd GPEngine::GetPerGPCSPLogLikelihoods(size_t start, size_t length) const {
  return log_likelihoods_.block(start, 0, length, log_likelihoods_.cols()) *
         site_pattern_weights_;
};

EigenVectorXd GPEngine::GetPerGPCSPComponentsOfFullLogMarginal() const {
  return GetPerGPCSPLogLikelihoods().array() +
         static_cast<double>(site_pattern_.SiteCount()) * q_.array().log();
};

EigenConstMatrixXdRef GPEngine::GetLogLikelihoodMatrix() const {
  return log_likelihoods_.block(0, 0, gpcsp_count_, log_likelihoods_.cols());
};

EigenConstVectorXdRef GPEngine::GetHybridMarginals() const {
  return hybrid_marginal_log_likelihoods_;
};

EigenConstVectorXdRef GPEngine::GetSBNParameters() const { return q_; };

EigenMatrixXd GPEngine::GetPLV(size_t plv_index) const { return plvs_.at(plv_index); }

void GPEngine::PrintPLV(size_t plv_idx) const {
  for (auto&& row : plvs_[plv_idx].rowwise()) {
    std::cout << row << std::endl;
  }
  std::cout << std::endl;
}

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

void GPEngine::InitializePLVsWithSitePatterns() {
  for (auto& plv : plvs_) {
    plv.setZero();
  }
  size_t taxon_idx = 0;
  for (const auto& pattern : site_pattern_.GetPatterns()) {
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative symbol!");
      if (symbol == MmappedNucleotidePLV::base_count_) {  // Gap character.
        plvs_.at(taxon_idx).col(site_idx).setConstant(1.);
      } else if (symbol < MmappedNucleotidePLV::base_count_) {
        plvs_.at(taxon_idx)(symbol, site_idx) = 1.;
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
  plvs_.at(plv_idx) /= pow(rescaling_threshold_, static_cast<double>(rescaling_count));
  rescaling_counts_(plv_idx) += rescaling_count;
}

void GPEngine::AssertPLVIsFinite(size_t plv_idx, const std::string& message) const {
  Assert(plvs_.at(plv_idx).array().isFinite().all(), message);
}

std::pair<double, double> GPEngine::PLVMinMax(size_t plv_idx) const {
  return {plvs_.at(plv_idx).minCoeff(), plvs_.at(plv_idx).maxCoeff()};
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

void GPEngine::BrentOptimization(const GPOperations::OptimizeBranchLength& op) {
  auto negative_log_likelihood = [this, &op](double log_branch_length) {
    SetTransitionMatrixToHaveBranchLength(exp(log_branch_length));
    PreparePerPatternLogLikelihoodsForGPCSP(op.rootward_, op.leafward_);
    return -per_pattern_log_likelihoods_.dot(site_pattern_weights_);
  };
  double current_log_branch_length = log(branch_lengths_(op.gpcsp_));
  double current_value = negative_log_likelihood(current_log_branch_length);
  const auto [log_branch_length, neg_log_likelihood] = Optimization::BrentMinimize(
      negative_log_likelihood, min_log_branch_length_, max_log_branch_length_,
      significant_digits_for_optimization_, max_iter_for_optimization_);

  // Numerical optimization sometimes yields new nllk > current nllk.
  // In this case, we reset the branch length to the previous value.
  if (neg_log_likelihood > current_value) {
    branch_lengths_(op.gpcsp_) = exp(current_log_branch_length);
  } else {
    branch_lengths_(op.gpcsp_) = exp(log_branch_length);
  }
}

void GPEngine::GradientAscentOptimization(
    const GPOperations::OptimizeBranchLength& op) {
  // See issue #278:
  Failwith(
      "Hasn't been tested after the change to log-space branch length optimization.");
  auto log_likelihood_and_derivative = [this, &op](double log_branch_length) {
    branch_lengths_(op.gpcsp_) = exp(log_branch_length);
    return this->LogLikelihoodAndDerivative(op);
  };
  const auto [log_branch_length, log_likelihood] = Optimization::GradientAscent(
      log_likelihood_and_derivative, log(branch_lengths_(op.gpcsp_)),
      relative_tolerance_for_optimization_, step_size_for_optimization_,
      min_log_branch_length_, max_iter_for_optimization_);
  std::ignore = log_likelihood;
  branch_lengths_(op.gpcsp_) = exp(log_branch_length);
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
      // Normalize the branch length total using the counts to get a mean branch length.
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
    quartet_root_plv_ = transition_matrix_ * plvs_.at(rootward_tip.plv_idx_);
    for (const auto& sister_tip : request.sister_tips_) {
      CheckRescaling(sister_tip.plv_idx_);
      // Form the PLV on the root side of the central edge.
      SetTransitionMatrixToHaveBranchLength(branch_lengths_[sister_tip.gpcsp_idx_]);
      quartet_r_s_plv_.array() =
          quartet_root_plv_.array() *
          (transition_matrix_ * plvs_.at(sister_tip.plv_idx_)).array();
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
            (transition_matrix_ * plvs_.at(rotated_tip.plv_idx_)).array();
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
               plvs_.at(sorted_tip.plv_idx_))
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
