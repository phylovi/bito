// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Enables hybrid marginalization on a three-taxon tree.

#ifndef SRC_TRIPOD_HYBRID_MARGINALIZER_HPP_
#define SRC_TRIPOD_HYBRID_MARGINALIZER_HPP_

#include "mmapped_plv.hpp"
#include "tripod_hybrid_request.hpp"

class TripodHybridMarginalizer {
 public:
  TripodHybridMarginalizer(const NucleotidePLVRefVector& plvs,
                           EigenConstVectorXdRef branch_lengths,
                           EigenVectorXd node_probabilities);

  //  double LogLikelihood(size_t plv_root_idx, size_t root_node_idx,
  //                       double root_branch_length);

  std::vector<double> Process(const TripodHybridRequest& request);

 private:
  const NucleotidePLVRefVector& plvs_;
  EigenConstVectorXdRef branch_lengths_;
  EigenVectorXd node_probabilities_;
  EigenMatrixXd tripod_root_plv_;
  EigenMatrixXd tripod_left_plv_;
  EigenMatrixXd tripod_right_plv_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TripodHybridMarginalizer") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TRIPOD_HYBRID_MARGINALIZER_HPP_
