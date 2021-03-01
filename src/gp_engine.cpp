// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"

#include "optimization.hpp"
#include "sugar.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t plv_count, size_t gpcsp_count,
                   const std::string& mmap_file_path, double rescaling_threshold,
                   EigenVectorXd sbn_prior,
                   EigenVectorXd unconditional_node_probabilities,
                   EigenVectorXd inverted_sbn_prior)
    : site_pattern_(std::move(site_pattern)),
      plv_count_(plv_count),
      rescaling_threshold_(rescaling_threshold),
      log_rescaling_threshold_(log(rescaling_threshold)),
      mmapped_master_plv_(mmap_file_path, plv_count_ * site_pattern_.PatternCount()),
      plvs_(mmapped_master_plv_.Subdivide(plv_count_)),
      q_(std::move(sbn_prior)),
      unconditional_node_probabilities_(std::move(unconditional_node_probabilities)),
      inverted_sbn_prior_(std::move(inverted_sbn_prior)) {
  Assert(plvs_.back().rows() == MmappedNucleotidePLV::base_count_ &&
             plvs_.back().cols() == site_pattern_.PatternCount(),
         "Didn't get the right shape of PLVs out of Subdivide.");
  rescaling_counts_.resize(plv_count_);
  rescaling_counts_.setZero();
  branch_lengths_.resize(gpcsp_count);
  branch_lengths_.setConstant(default_branch_length_);
  log_marginal_likelihood_.resize(site_pattern_.PatternCount());
  log_marginal_likelihood_.setConstant(DOUBLE_NEG_INF);
  log_likelihoods_.resize(gpcsp_count, site_pattern_.PatternCount());

  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

  quartet_root_plv_ = plvs_.at(0);
  quartet_root_plv_.setZero();
  quartet_r_s_plv_ = quartet_root_plv_;
  quartet_q_s_plv_ = quartet_root_plv_;
  quartet_r_sorted_plv_ = quartet_root_plv_;

  InitializePLVsWithSitePatterns();
}

void GPEngine::operator()(const GPOperations::ZeroPLV& op) {
  plvs_.at(op.dest_).setZero();
  rescaling_counts_(op.dest_) = 0;
}

void GPEngine::operator()(const GPOperations::SetToStationaryDistribution& op) {
  auto& plv = plvs_.at(op.dest_);
  for (size_t row_idx = 0; row_idx < plv.rows(); ++row_idx) {
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

EigenVectorXd NormalizedPosteriorOfUnnormalized(EigenVectorXd unnormalized_posterior) {
  const double log_norm = NumericalUtils::LogSum(unnormalized_posterior);
  unnormalized_posterior.array() -= log_norm;
  return unnormalized_posterior.array().exp();
}

void GPEngine::operator()(const GPOperations::UpdateSBNProbabilities& op) {
  const size_t range_length = op.stop_ - op.start_;
  if (range_length == 1) {
    q_(op.start_) = 1.;
  } else {
    EigenVectorXd log_likelihoods = GetPerGPCSPLogLikelihoods(op.start_, range_length);
    EigenVectorXd log_prior = q_.segment(op.start_, range_length).array().log();
    q_.segment(op.start_, range_length) =
        NormalizedPosteriorOfUnnormalized(log_likelihoods + log_prior);
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
  Assert(branch_lengths_.size() == branch_lengths.size(),
         "Size mismatch in GPEngine::SetBranchLengths.");
  branch_lengths_ = std::move(branch_lengths);
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

EigenVectorXd GPEngine::GetBranchLengths() const { return branch_lengths_; };

EigenVectorXd GPEngine::GetPerGPCSPLogLikelihoods() const {
  return log_likelihoods_ * site_pattern_weights_;
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
  return log_likelihoods_;
};

EigenConstVectorXdRef GPEngine::GetSBNParameters() const { return q_; };

void GPEngine::PrintPLV(size_t plv_idx) {
  for (const auto& row : plvs_[plv_idx].rowwise()) {
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
  branch_lengths_(op.gpcsp_) = exp(log_branch_length);
}

void GPEngine::HotStartBranchLengths(const RootedTreeCollection& tree_collection,
                                     const BitsetSizeMap& indexer) {
  const auto leaf_count = tree_collection.TaxonCount();
  const size_t default_index = branch_lengths_.size();
  branch_lengths_.setZero();
  Eigen::VectorXi gpcsp_counts = Eigen::VectorXi::Zero(branch_lengths_.size());
  // Set the branch length vector to be the total of the branch lengths for each PCSP,
  // and count the number of times we have seen each PCSP (into gpcsp_counts).
  for (const auto& tree : tree_collection.Trees()) {
    tree.Topology()->RootedPCSPPreorder(
        [&leaf_count, &default_index, &indexer, &tree, &gpcsp_counts, this](
            const Node* sister_node, const Node* focal_node, const Node* child0_node,
            const Node* child1_node) {
          Bitset gpcsp_bitset =
              SBNMaps::PCSPBitsetOf(leaf_count, sister_node, false, focal_node, false,
                                    child0_node, false, child1_node, false);
          const auto gpcsp_idx = AtWithDefault(indexer, gpcsp_bitset, default_index);
          if (gpcsp_idx != default_index) {
            branch_lengths_(gpcsp_idx) += tree.BranchLength(focal_node);
            gpcsp_counts(gpcsp_idx)++;
          }
        });
  }
  for (Eigen::Index gpcsp_idx = 0; gpcsp_idx < gpcsp_counts.size(); ++gpcsp_idx) {
    if (gpcsp_counts(gpcsp_idx) == 0) {
      branch_lengths_(gpcsp_idx) = default_branch_length_;
    } else {
      // Normalize the branch length total using the counts to get a mean branch length.
      branch_lengths_(gpcsp_idx) /= static_cast<double>(gpcsp_counts(gpcsp_idx));
    }
  }
}

// #323 rescaling factor
// #323 transpose for down the tree
std::vector<double> GPEngine::ProcessQuartetHybridRequest(
    const QuartetHybridRequest& request) {
  std::vector<double> result;
  for (const auto& rootward_tip : request.rootward_tips_) {
    const double rootward_tip_prior =
        unconditional_node_probabilities_[rootward_tip.tip_node_id_];
    const double log_rootward_tip_prior = log(rootward_tip_prior);
    SetTransitionMatrixToHaveBranchLength(branch_lengths_[rootward_tip.gpcsp_idx_]);
    quartet_root_plv_ = transition_matrix_ * plvs_.at(rootward_tip.plv_idx_);
    for (const auto& sister_tip : request.sister_tips_) {
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
        // Form the PLV on the root side of the sorted edge.
        SetTransitionMatrixToHaveBranchLength(branch_lengths_[rotated_tip.gpcsp_idx_]);
        quartet_r_sorted_plv_.array() =
            quartet_q_s_plv_.array() *
            (transition_matrix_ * plvs_.at(rotated_tip.plv_idx_)).array();
        for (const auto& sorted_tip : request.sorted_tips_) {
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
  return result;
}
