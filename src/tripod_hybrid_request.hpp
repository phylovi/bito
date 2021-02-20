// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Enables hybrid marginalization on a three-taxon tree.

#ifndef SRC_TRIPOD_HYBRID_REQUEST_HPP_
#define SRC_TRIPOD_HYBRID_REQUEST_HPP_

#include <vector>

struct PLVPCSPPair {
  size_t plv_id;
  size_t pcsp_id;
};

struct TripodHybridRequest {
  std::vector<PLVPCSPPair> rootward_pairs;
  std::vector<PLVPCSPPair> rotated_pairs;
  std::vector<PLVPCSPPair> sorted_pairs;

  size_t central_pcsp;
};

#endif  // SRC_TRIPOD_HYBRID_REQUEST_HPP_
