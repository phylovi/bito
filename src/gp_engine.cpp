// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"
#include "optimization.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t plv_count,
                   size_t continous_parameter_count, std::string mmap_file_path)
    : site_pattern_(std::move(site_pattern)),
      plv_count_(plv_count),
      mmapped_master_plv_(mmap_file_path, plv_count_ * site_pattern_.PatternCount()) {
  Assert(plv_count_ > 0, "Zero PLV count in constructor of GPEngine.");
  plvs_ = mmapped_master_plv_.Subdivide(plv_count_);
  Assert(plvs_.size() == plv_count_,
         "Didn't get the right number of PLVs out of Subdivide.");
  Assert(plvs_.back().rows() == MmappedNucleotidePLV::base_count_ &&
             plvs_.back().cols() == site_pattern_.PatternCount(),
         "Didn't get the right shape of PLVs out of Subdivide.");
  branch_lengths_.resize(continous_parameter_count);
  branch_lengths_.setOnes();
  log_likelihoods_.resize(continous_parameter_count);
  q_.resize(continous_parameter_count);

  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

  InitializePLVsWithSitePatterns();
}

void GPEngine::operator()(const GPOperations::Zero& op) {
  plvs_.at(op.dest_idx).setZero();
}

void GPEngine::operator()(const GPOperations::SetToStationaryDistribution& op) {
  auto& plv = plvs_.at(op.dest_idx);
  for (size_t row_idx = 0; row_idx < plv.rows(); ++row_idx) {
    plv.row(row_idx).array() = stationary_distribution_(row_idx);
  }
}

void GPEngine::operator()(const GPOperations::EvolvePLVWeightedBySBNParameter& op) {
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(op.gpcsp_idx));
  plvs_.at(op.dest_idx) += q_(op.gpcsp_idx) * transition_matrix_ * plvs_.at(op.src_idx);
}

void GPEngine::operator()(const GPOperations::IncrementMarginalLikelihood& op) {
  per_pattern_log_likelihoods_ =
      (plvs_.at(op.stationary_idx).transpose() * plvs_.at(op.p_idx))
          .diagonal()
          .array()
          .log();
  log_likelihoods_[op.gpcsp_idx] =
      log(q_(op.p_idx)) + per_pattern_log_likelihoods_.dot(site_pattern_weights_);
  log_marginal_likelihood_ =
      NumericalUtils::LogAdd(log_marginal_likelihood_, log_likelihoods_[op.gpcsp_idx]);
}

void GPEngine::operator()(const GPOperations::Multiply& op) {
  plvs_.at(op.dest_idx).array() =
      plvs_.at(op.src1_idx).array() * plvs_.at(op.src2_idx).array();
}

void GPEngine::operator()(const GPOperations::Likelihood& op) {
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(op.dest_idx));
  auto result = plvs_.at(op.parent_idx).transpose() *
                (transition_matrix_ * plvs_.at(op.child_idx));
  per_pattern_log_likelihoods_ = result.diagonal().array().log();
  log_likelihoods_[op.dest_idx] =
      log(q_[op.dest_idx]) + per_pattern_log_likelihoods_.dot(site_pattern_weights_);
}

void GPEngine::operator()(const GPOperations::OptimizeBranchLength& op) {
  GradientAscentOptimization(op);
}

void GPEngine::operator()(const GPOperations::UpdateSBNProbabilities& op) {
  size_t range_length = op.stop_idx - op.start_idx;
  if (range_length == 1) {
    q_(op.start_idx) = 1.;
    return;
  }

  auto segment = log_likelihoods_.segment(op.start_idx, op.stop_idx - op.start_idx);
  double log_norm = NumericalUtils::LogSum(segment);
  segment.array() -= log_norm;
  q_.segment(op.start_idx, op.stop_idx - op.start_idx) = segment.array().exp();
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
  Eigen::Vector4d diagonal_vector = (branch_length * eigenvalues_).array().exp();
  diagonal_matrix_.diagonal() = diagonal_vector;
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
  diagonal_matrix_.diagonal() = eigenvalues_.array() * diagonal_vector.array();
  derivative_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void GPEngine::SetTransitionMatrixToHaveBranchLengthAndTranspose(double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ =
      inverse_eigenmatrix_.transpose() * diagonal_matrix_ * eigenmatrix_.transpose();
}

void GPEngine::PrintPLV(size_t plv_idx) {
  for (auto row : plvs_[plv_idx].rowwise()) {
    std::cout << row << std::endl;
  }
  std::cout << std::endl;
}

DoublePair GPEngine::LogLikelihoodAndDerivative(
    const GPOperations::OptimizeBranchLength& op) {
  SetTransitionAndDerivativeMatricesToHaveBranchLength(branch_lengths_(op.gpcsp_idx));

  // The per-site likelihood derivative is calculated in the same way as the per-site
  // likelihood, but using the derivative matrix instead of the transition matrix.
  PreparePerPatternLikelihoodDerivatives(op.rootward_idx, op.leafward_idx);
  PreparePerPatternLikelihoods(op.rootward_idx, op.leafward_idx);

  per_pattern_log_likelihoods_ = per_pattern_likelihoods_.array().log();
  // The phylogenetic component of the likelihood is weighted with the number of times
  // we see the site patterns.
  // TODO Seong-- I'm trying to understand the role of q_ here. Is this as a placeholder
  // for the prior?
  const double log_likelihood =
      log(q_(op.gpcsp_idx)) + per_pattern_log_likelihoods_.dot(site_pattern_weights_);

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

void GPEngine::BrentOptimization(const GPOperations::OptimizeBranchLength& op) {
  auto negative_log_likelihood = [this, &op](double branch_length) {
    SetTransitionMatrixToHaveBranchLength(branch_length);
    auto result = plvs_.at(op.rootward_idx).transpose() *
                  (transition_matrix_ * plvs_.at(op.leafward_idx));
    per_pattern_log_likelihoods_ = result.diagonal().array().log();
    return -(log(q_[op.gpcsp_idx]) +
             per_pattern_log_likelihoods_.dot(site_pattern_weights_));
  };
  double current_branch_length = branch_lengths_(op.gpcsp_idx);
  double current_value = negative_log_likelihood(current_branch_length);
  auto [branch_length, neg_log_likelihood] = Optimization::BrentMinimize(
      negative_log_likelihood, min_branch_length_, max_branch_length_,
      significant_digits_for_optimization_, max_iter_for_optimization_);

  // Numerical optimization sometimes yields new nllk > current nllk.
  // In this case, we reset the branch length to the previous value.
  if (neg_log_likelihood > current_value) {
    branch_lengths_(op.gpcsp_idx) = current_branch_length;
  } else {
    branch_lengths_(op.gpcsp_idx) = branch_length;
  }
}

void GPEngine::GradientAscentOptimization(
    const GPOperations::OptimizeBranchLength& op) {
  auto log_likelihood_and_derivative = [this, &op](double branch_length) {
    branch_lengths_(op.gpcsp_idx) = branch_length;
    return this->LogLikelihoodAndDerivative(op);
  };
  auto [branch_length, log_likelihood] = Optimization::GradientAscent(
      log_likelihood_and_derivative, branch_lengths_(op.gpcsp_idx),
      relative_tolerance_for_optimization_, step_size_for_optimization_,
      min_branch_length_, max_iter_for_optimization_);
  branch_lengths_(op.gpcsp_idx) = branch_length;
}
