// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"
#include "optimization.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t gpcsp_count)
    : site_pattern_(std::move(site_pattern)) {
  auto plv_count = site_pattern_.PatternCount() + gpcsp_count;
  plvs_ = std::vector<NucleotidePLV>(
      plv_count, NucleotidePLV::Zero(4, site_pattern_.PatternCount()));
  branch_lengths_.resize(gpcsp_count);
  log_likelihoods_.resize(gpcsp_count);
  q_.resize(gpcsp_count);

  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

  InitializePLVsWithSitePatterns();
}

void GPEngine::InitializePLVsWithSitePatterns() {
  size_t taxon_idx = 0;
  for (const auto& pattern : site_pattern_.GetPatterns()) {
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative symbol!");
      if (symbol == 4) {  // Gap character.
        plvs_.at(taxon_idx).col(site_idx).setConstant(1.);
      } else if (symbol < 4) {
        plvs_.at(taxon_idx)(symbol, site_idx) = 1.;
      }
      site_idx++;
    }
    taxon_idx++;
  }
}

void GPEngine::operator()(const GPOperations::Zero& op) {
  plvs_.at(op.dest_idx).setZero();
}

void GPEngine::operator()(const GPOperations::SetToStationaryDistribution& op) {
  auto& plv = plvs_.at(op.dest_idx);
  for (size_t row_idx = 0; row_idx < 4; ++row_idx) {
    plv.row(row_idx).array() = stationary_distribution_(row_idx);
  }
}

void GPEngine::operator()(const GPOperations::WeightedSumAccumulate& op) {
  plvs_.at(op.dest_idx) += q_(op.q_idx) * plvs_.at(op.src_idx);
}

void GPEngine::operator()(const GPOperations::Multiply& op) {
  plvs_.at(op.dest_idx).array() =
      plvs_.at(op.src1_idx).array() * plvs_.at(op.src2_idx).array();
}

void GPEngine::operator()(const GPOperations::Likelihood& op) {
  log_likelihoods_(op.dest_idx) = LogLikelihood(op.src1_idx, op.src2_idx);
}

void GPEngine::operator()(const GPOperations::EvolveRootward& op) {
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(op.branch_length_idx));
  plvs_.at(op.dest_idx) = transition_matrix_ * plvs_.at(op.src_idx);
}

void GPEngine::operator()(const GPOperations::EvolveLeafward& op) {
  SetTransitionMatrixToHaveBranchLengthAndTranspose(
      branch_lengths_(op.branch_length_idx));
  plvs_.at(op.dest_idx) = transition_matrix_ * plvs_.at(op.src_idx);
}

DoublePair GPEngine::LogLikelihoodAndDerivative(
    const GPOperations::OptimizeRootward& op) {
  SetTransitionAndDerivativeMatricesToHaveBranchLength(
      branch_lengths_(op.branch_length_idx));
  // The per-site likelihood derivative is calculated in the same way as the per-site
  // likelihood, but using the derivative matrix instead of the transition matrix.
  plvs_.at(op.dest_idx) = derivative_matrix_ * plvs_.at(op.leafward_idx);
  PreparePerPatternLikelihoodDerivatives(op.rootward_idx, op.dest_idx);
  plvs_.at(op.dest_idx) = transition_matrix_ * plvs_.at(op.leafward_idx);
  PreparePerPatternLikelihoods(op.rootward_idx, op.dest_idx);
  return LogLikelihoodAndDerivativeFromPreparations();
}

void GPEngine::BrentOptimization(const GPOperations::OptimizeRootward& op) {
  auto negative_log_likelihood = [this, &op](double branch_length) {
    SetTransitionMatrixToHaveBranchLength(branch_length);
    plvs_.at(op.dest_idx) = transition_matrix_ * plvs_.at(op.leafward_idx);
    return -LogLikelihood(op.rootward_idx, op.dest_idx);
  };
  auto [branch_length, neg_log_likelihood] = Optimization::BrentMinimize(
      negative_log_likelihood, branch_length_min_, branch_length_max_,
      significant_digits_for_optimization_, max_iter_for_optimization_);
  branch_lengths_(op.branch_length_idx) = branch_length;
  log_likelihoods_(op.branch_length_idx) = -neg_log_likelihood;
}

DoublePair GradientAscent(std::function<DoublePair(double)> f_and_f_prime, double x,
                          const double tolerance, const double step_size,
                          const double min_x, const size_t max_iter) {
  size_t iter_idx = 0;
  // std::cout << "x, f_x, f_prime_x\n";
  while (true) {
    auto [f_x, f_prime_x] = f_and_f_prime(x);
    // std::cout << "gradient ascent " << iter_idx << " " << std::vector{x, f_x,
    // f_prime_x}
    //           << std::endl;
    const double new_x = x + f_prime_x * step_size;
    x = std::max(new_x, min_x);
    if (fabs(f_prime_x) < fabs(f_x) * tolerance || iter_idx >= max_iter) {
      return {x, f_x};
    }
    ++iter_idx;
  }
}

void GPEngine::GradientAscentOptimization(const GPOperations::OptimizeRootward& op) {
  auto log_likelihood_and_derivative = [this, &op](double branch_length) {
    branch_lengths_(op.branch_length_idx) = branch_length;
    return this->LogLikelihoodAndDerivative(op);
  };
  auto [branch_length, log_likelihood] = GradientAscent(
      log_likelihood_and_derivative, branch_lengths_(op.branch_length_idx),
      relative_tolerance_for_optimization_, step_size_for_optimization_,
      branch_length_min_, max_iter_for_optimization_);
  branch_lengths_(op.branch_length_idx) = branch_length;
  log_likelihoods_(op.branch_length_idx) = log_likelihood;
}

void GPEngine::operator()(const GPOperations::OptimizeRootward& op) {
  auto starting_branch_length = branch_lengths_(op.branch_length_idx);
  std::cout << "starting branch length: " << starting_branch_length << std::endl;
  GradientAscentOptimization(op);
  std::cout << "after gradient ascent: " << branch_lengths_(op.branch_length_idx)
            << std::endl;
  branch_lengths_(op.branch_length_idx) = starting_branch_length;
  BrentOptimization(op);
  std::cout << "after brent: " << branch_lengths_(op.branch_length_idx) << std::endl;
}

void GPEngine::operator()(const GPOperations::OptimizeLeafward& op) {
  Failwith("OptimizeRootward unimplemented for now.");
}

void GPEngine::operator()(const GPOperations::UpdateSBNProbabilities& op) {
  Failwith("UpdateSBNProbabilities unimplemented for now.");
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
  diagonal_matrix_.diagonal() = diagonal_vector.array() * eigenvalues_.array();
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
