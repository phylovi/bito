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

  void operator()(const GPOperations::Zero& operation) { std::cout << "zero process"; }
  void operator()(const GPOperations::WeightedSumAccumulate& operation) {}
  void operator()(const GPOperations::Multiply& operation) {}
  void operator()(const GPOperations::Likelihood& operation) {}
  void operator()(const GPOperations::EvolveRootward& operation) {}
  void operator()(const GPOperations::EvolveLeafward& operation) {}
  void operator()(const GPOperations::OptimizeRootward& operation) {}
  void operator()(const GPOperations::OptimizeLeafward& operation) {}
  void operator()(const GPOperations::UpdateSBNProbabilities& operation) {}

  void ProcessOperations(GPOperationVector operations) {
    for (const auto& operation : operations) {
      std::visit(*this, operation);
    }
  }

 private:
  SitePattern site_pattern_;
  size_t pcss_count_;
  std::vector<NucleotidePLV> plvs_;
};

#endif  // SRC_GP_ENGINE_HPP_
