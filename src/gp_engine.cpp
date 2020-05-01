// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"
#include "optimization.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t pcss_count)
    : site_pattern_(std::move(site_pattern)), pcss_count_(pcss_count) {
  auto plv_count = site_pattern_.PatternCount() + pcss_count_;
  plvs_ = std::vector<NucleotidePLV>(
      plv_count, NucleotidePLV::Zero(4, site_pattern_.PatternCount()));
  branch_lengths_.resize(pcss_count_);
  log_likelihoods_.resize(pcss_count_);
  q_.resize(pcss_count_);

  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

  InitializePLVsWithSitePatterns();
}

void GPEngine::InitializePLVsWithSitePatterns() {
  size_t taxon_idx = 0;
  for (const auto &pattern : site_pattern_.GetPatterns()) {
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

void GPEngine::operator()(const GPOperations::OptimizeRootward& op) {
  auto to_optimize = [this, &op](double branch_length) {
    SetTransitionMatrixToHaveBranchLength(branch_length);
    plvs_.at(op.dest_idx) = transition_matrix_ * plvs_.at(op.leafward_idx);
    return LogLikelihood(op.rootward_idx, op.dest_idx);
  };
  // TODO can we start optimization somewhere smart?
  auto [branch_length, log_likelihood] =
      Optimization::BrentMinimize(to_optimize, 0, 5, 4, 20);
  branch_lengths_[op.branch_length_idx] = branch_length;
  log_likelihoods_[op.dest_idx] = log_likelihood;
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
