// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Stores a "request" for a tripod hybrid marginal calculation.

#ifndef SRC_TRIPOD_HYBRID_REQUEST_HPP_
#define SRC_TRIPOD_HYBRID_REQUEST_HPP_

#include <iostream>
#include <vector>

struct TripodTip {
  constexpr TripodTip(size_t tip_node_id, size_t plv_idx, size_t gpcsp_idx)
      : tip_node_id_(tip_node_id), plv_idx_(plv_idx), gpcsp_idx_(gpcsp_idx){};

  size_t tip_node_id_;
  size_t plv_idx_;
  size_t gpcsp_idx_;
};

using TripodTipVector = std::vector<TripodTip>;

struct TripodHybridRequest {
  TripodHybridRequest(size_t central_gpcsp_idx, TripodTipVector rootward_pairs,
                      TripodTipVector rotated_pairs, TripodTipVector sorted_pairs)
      : central_gpcsp_idx_(central_gpcsp_idx),
        rootward_pairs_(std::move(rootward_pairs)),
        rotated_pairs_(std::move(rotated_pairs)),
        sorted_pairs_(std::move(sorted_pairs)){};

  size_t central_gpcsp_idx_;
  TripodTipVector rootward_pairs_;
  TripodTipVector rotated_pairs_;
  TripodTipVector sorted_pairs_;
};

std::ostream& operator<<(std::ostream& os, TripodTip const& plv_pcsp);
std::ostream& operator<<(std::ostream& os, TripodTipVector const& plv_pcsp_vector);
std::ostream& operator<<(std::ostream& os, TripodHybridRequest const& request);

#endif  // SRC_TRIPOD_HYBRID_REQUEST_HPP_
