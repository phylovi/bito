// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A crappy version of Engine that only does JC, but does do our GP calculations.

#ifndef SRC_GP_ENGINE_HPP_
#define SRC_GP_ENGINE_HPP_

#include "eigen_sugar.hpp"
#include "gp_operation.hpp"
#include "site_pattern.hpp"

using NucleotidePLV = Eigen::Matrix<double, Eigen::Dynamic, 4, Eigen::RowMajor>;

class GPEngine {
 public:
  GPEngine(SitePattern site_pattern, size_t pcss_count);

  void operator()(const GPOperations::Zero& op) { plvs_[op.dest_idx].setZero(); }
  void operator()(const GPOperations::WeightedSumAccumulate& op) {
    plvs_[op.dest_idx] += q_[op.q_idx] * plvs_[op.src_idx];
  }
  void operator()(const GPOperations::Multiply& op) {
    plvs_[op.dest_idx].array() =
        plvs_[op.src1_idx].array() * plvs_[op.src2_idx].array();
  }
  void operator()(const GPOperations::Likelihood& op) {}
  void operator()(const GPOperations::EvolveRootward& op) {}
  void operator()(const GPOperations::EvolveLeafward& op) {}
  void operator()(const GPOperations::OptimizeRootward& op) {}
  void operator()(const GPOperations::OptimizeLeafward& op) {}
  void operator()(const GPOperations::UpdateSBNProbabilities& op) {}

  void ProcessOperations(GPOperationVector operations) {
    for (const auto& operation : operations) {
      std::visit(*this, operation);
    }
  }

 private:
  SitePattern site_pattern_;
  size_t pcss_count_;
  std::vector<NucleotidePLV> plvs_;
  EigenVectorXd branch_lengths_;
  EigenVectorXd likelihoods_;
  EigenVectorXd q_;

  void InitializePLVsWithSitePatterns();
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("GPEngine") {}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_GP_ENGINE_HPP_
