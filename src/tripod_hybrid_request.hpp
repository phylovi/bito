// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Stores a "request" for a tripod hybrid marginal calculation.

#ifndef SRC_TRIPOD_HYBRID_REQUEST_HPP_
#define SRC_TRIPOD_HYBRID_REQUEST_HPP_

#include <iostream>
#include <vector>

struct PLVPCSPPair {
  constexpr PLVPCSPPair(size_t plv_idx, size_t gpcsp_idx)
      : plv_idx_(plv_idx), gpcsp_idx_(gpcsp_idx){};

  size_t plv_idx_;
  size_t gpcsp_idx_;
};

using PLVPCSPPairVector = std::vector<PLVPCSPPair>;

struct TripodHybridRequest {
  TripodHybridRequest(size_t central_gpcsp_idx, PLVPCSPPairVector rootward_pairs,
                      PLVPCSPPairVector rotated_pairs, PLVPCSPPairVector sorted_pairs)
      : central_gpcsp_idx_(central_gpcsp_idx),
        rootward_pairs_(std::move(rootward_pairs)),
        rotated_pairs_(std::move(rotated_pairs)),
        sorted_pairs_(std::move(sorted_pairs)){};

  size_t central_gpcsp_idx_;
  std::vector<PLVPCSPPair> rootward_pairs_;
  std::vector<PLVPCSPPair> rotated_pairs_;
  std::vector<PLVPCSPPair> sorted_pairs_;
};

std::ostream& operator<<(std::ostream& os, PLVPCSPPair const& plv_pcsp);
std::ostream& operator<<(std::ostream& os, PLVPCSPPairVector const& plv_pcsp_vector);
std::ostream& operator<<(std::ostream& os, TripodHybridRequest const& request);

#endif  // SRC_TRIPOD_HYBRID_REQUEST_HPP_
