// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t pcss_count)
    : site_pattern_(std::move(site_pattern)), pcss_count_(pcss_count) {
  auto plv_count = site_pattern_.PatternCount() + pcss_count_;
  plvs_ = std::vector<NucleotidePLV>(
      plv_count, NucleotidePLV::Zero(4, site_pattern_.PatternCount()));
  branch_lengths_.resize(pcss_count_);
  likelihoods_.resize(pcss_count_);
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

void GPEngine::ProcessOperations(GPOperationVector operations) {
  for (const auto& operation : operations) {
    std::visit(*this, operation);
  }
}

void GPEngine::SetBranchLengthForTransitionMatrix(double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void GPEngine::PrintPLV(size_t plv_idx) {
  for (auto row : plvs_[plv_idx].rowwise()) {
    std::cout << row << std::endl;
  }
  std::cout << std::endl;
}
