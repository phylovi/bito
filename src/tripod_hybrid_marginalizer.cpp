// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tripod_hybrid_marginalizer.hpp"

TripodHybridMarginalizer::TripodHybridMarginalizer(const NucleotidePLVRefVector& plvs,
                                                   EigenConstVectorXdRef branch_lengths)
    : plvs_(plvs), branch_lengths_(branch_lengths) {}
