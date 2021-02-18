// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Enables hybrid marginalization on a three-taxon tree.

#ifndef SRC_TRIPOD_HYBRID_MARGINALIZER_HPP_
#define SRC_TRIPOD_HYBRID_MARGINALIZER_HPP_

#include "mmapped_plv.hpp"

class TripodHybridMarginalizer {
  TripodHybridMarginalizer(const NucleotidePLVRefVector& plvs_,
                           EigenConstVectorXdRef branch_lengths_);

 private:
  const NucleotidePLVRefVector& plvs_;
  EigenConstVectorXdRef branch_lengths_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TripodHybridMarginalizer") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TRIPOD_HYBRID_MARGINALIZER_HPP_
