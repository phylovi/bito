// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"

#include "optimization.hpp"
#include "sugar.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t plv_count, size_t gpcsp_count,
                   const std::string& mmap_file_path, double rescaling_threshold)
    : site_pattern_(std::move(site_pattern)),
      plv_count_(plv_count),
      rescaling_threshold_(rescaling_threshold),
      log_rescaling_threshold_(log(rescaling_threshold)),
      mmapped_master_plv_(mmap_file_path, plv_count_ * site_pattern_.PatternCount()) {
  Assert(plv_count_ > 0, "Zero PLV count in constructor of GPEngine.");
  plvs_ = mmapped_master_plv_.Subdivide(plv_count_);
  Assert(plvs_.size() == plv_count_,
         "Didn't get the right number of PLVs out of Subdivide.");
  Assert(plvs_.back().rows() == MmappedNucleotidePLV::base_count_ &&
             plvs_.back().cols() == site_pattern_.PatternCount(),
         "Didn't get the right shape of PLVs out of Subdivide.");
  rescaling_counts_ = EigenVectorXi::Zero(plv_count_);
  branch_lengths_.resize(gpcsp_count);
  branch_lengths_.setConstant(default_branch_length_);
  log_marginal_likelihood_.resize(site_pattern_.PatternCount());
  log_marginal_likelihood_.setConstant(DOUBLE_NEG_INF);
  log_likelihoods_.resize(gpcsp_count, site_pattern_.PatternCount());
  q_.resize(gpcsp_count);

  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

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

void GPEngine::operator()(const GPOperations::IncrementMarginalLikelihood& op) {
  Assert(rescaling_counts_(op.stationary_) == 0,
         "Surprise! Rescaled stationary distribution in IncrementMarginalLikelihood");
  log_likelihoods_.row(op.rootsplit_) =
      (plvs_.at(op.stationary_).transpose() * plvs_.at(op.p_))
          .diagonal()
          .array()
          .log() +
      LogRescalingFor(op.p_);

  // Perform LogAdd per site.
  log_marginal_likelihood_ = NumericalUtils::LogAddVectors(
      log_marginal_likelihood_, log_likelihoods_.row(op.rootsplit_));
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

void GPEngine::operator()(const GPOperations::UpdateSBNProbabilities& op) {
  const size_t range_length = op.stop_ - op.start_;
  if (range_length == 1) {
    q_(op.start_) = 1.;
    return;
  }
  // else
  EigenVectorXd log_likelihoods = GetPerGPCSPLogLikelihoods(op.start_, range_length);
  EigenVectorXd log_prior = q_.segment(op.start_, range_length).array().log();
  EigenVectorXd unnormalized_posterior = log_likelihoods + log_prior;
  const double log_norm = NumericalUtils::LogSum(unnormalized_posterior);
  unnormalized_posterior.array() -= log_norm;
  q_.segment(op.start_, range_length) = unnormalized_posterior.array().exp();
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

void GPEngine::PrintPLV(size_t plv_idx) {
  for (const auto& row : plvs_[plv_idx].rowwise()) {
    std::cout << row << std::endl;
  }
  std::cout << std::endl;
}

EigenVectorXd GPEngine::GetPerGPCSPLogLikelihoods(size_t start, size_t length) const {
  return log_likelihoods_.block(start, 0, length, log_likelihoods_.cols()) *
         site_pattern_weights_;
};

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
  // TODO why are we getting radically different marginals for different edges? What if
  // we turn off branch optimization?
  std::cout << "\nNEW BRANCH\ntrue log marginal: "
            << (log_marginal_likelihood_.array() * site_pattern_weights_.array()).sum()
            << std::endl;
  Assert(rescaling_counts_.maxCoeff() == 0, "We are rescaling in BrentOptimization.");
  PrepareUnrescaledPerPatternLikelihoods(op.rootward_, op.leafward_);
  // 5: C := F - D
  per_pattern_marginal_constant_ = log_marginal_likelihood_.array().exp() -
                                   per_pattern_likelihoods_.array() * q_(op.gpcsp_);
  std::cout << "log marg: " << log_marginal_likelihood_.array().exp() << std::endl;
  std::cout << "like: " << per_pattern_likelihoods_ << std::endl;
  std::cout << "const: " << per_pattern_marginal_constant_ << std::endl;
  std::cout << "q: " << q_(op.gpcsp_) << std::endl;
  // What we optimize.
  auto negative_log_likelihood = [this, &op](double log_branch_length) {
    SetTransitionMatrixToHaveBranchLength(exp(log_branch_length));
    PrepareUnrescaledPerPatternLikelihoods(op.rootward_, op.leafward_);
    // 6: f is log(D + C)
    per_pattern_log_likelihoods_ =
        (per_pattern_likelihoods_ * q_(op.gpcsp_) + per_pattern_marginal_constant_)
            .array()
            .log();
    std::cout << "branch_length: " << exp(log_branch_length) << std::endl;
    std::cout << "result: " << -per_pattern_log_likelihoods_.dot(site_pattern_weights_)
              << std::endl;
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
    // // 9: F = D + C
    // log_marginal_likelihood_ = (per_pattern_likelihoods_.array() * q_(op.gpcsp_) +
    //                             per_pattern_marginal_constant_.array())
    //                                .log();
  }
}

void GPEngine::PartialBrentOptimization(const GPOperations::OptimizeBranchLength& op) {
  auto negative_log_likelihood = [this, &op](double log_branch_length) {
    SetTransitionMatrixToHaveBranchLength(exp(log_branch_length));
    PreparePerPatternLogLikelihoodsForGPCSP(op.rootward_, op.leafward_);
    std::cout << "branch_length: " << exp(log_branch_length) << std::endl;
    std::cout << "result: " << -per_pattern_log_likelihoods_.dot(site_pattern_weights_)
              << std::endl;
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
