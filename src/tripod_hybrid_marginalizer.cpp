// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tripod_hybrid_marginalizer.hpp"

TripodHybridMarginalizer::TripodHybridMarginalizer(const NucleotidePLVRefVector& plvs,
                                                   EigenConstVectorXdRef branch_lengths,
                                                   EigenVectorXd node_probabilities)
    : plvs_(plvs),
      branch_lengths_(branch_lengths),
      node_probabilities_(std::move(node_probabilities)) {
  Assert(!plvs_.empty(),
         "Need at least one PLV in constructor of TripodHybridMarginalizer.");
  root_plv_ = plvs_.at(0);
  root_plv_.setZero();
  left_plv_ = root_plv_;
  right_plv_ = root_plv_;
}

std::vector<double> TripodHybridMarginalizer::Process(
    const TripodHybridRequest& request) {
  root_plv_(0, 0) = 1.;
  std::cout << root_plv_.row(0) << std::endl;
  std::cout << left_plv_.row(0) << std::endl;
  return {};
}
