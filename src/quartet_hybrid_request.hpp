// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Stores a "request" for a quartet hybrid marginal calculation.

#ifndef SRC_QUARTET_HYBRID_REQUEST_HPP_
#define SRC_QUARTET_HYBRID_REQUEST_HPP_

#include <iostream>
#include <vector>

struct QuartetTip {
  constexpr QuartetTip(size_t tip_node_id, size_t plv_idx, size_t gpcsp_idx)
      : tip_node_id_(tip_node_id), plv_idx_(plv_idx), gpcsp_idx_(gpcsp_idx){};
  size_t tip_node_id_;
  size_t plv_idx_;
  size_t gpcsp_idx_;
};

using QuartetTipVector = std::vector<QuartetTip>;

struct QuartetHybridRequest {
  QuartetHybridRequest(size_t central_gpcsp_idx, QuartetTipVector rootward_tips,
                       QuartetTipVector sister_tips, QuartetTipVector rotated_tips,
                       QuartetTipVector sorted_tips)
      : central_gpcsp_idx_(central_gpcsp_idx),
        rootward_tips_(std::move(rootward_tips)),
        sister_tips_(std::move(sister_tips)),
        rotated_tips_(std::move(rotated_tips)),
        sorted_tips_(std::move(sorted_tips)){};

  // Are each of the tip vectors non-empty?
  bool IsFullyFormed() const;

  size_t central_gpcsp_idx_;
  QuartetTipVector rootward_tips_;
  QuartetTipVector sister_tips_;
  QuartetTipVector rotated_tips_;
  QuartetTipVector sorted_tips_;
};

std::ostream& operator<<(std::ostream& os, QuartetTip const& plv_pcsp);
std::ostream& operator<<(std::ostream& os, QuartetTipVector const& plv_pcsp_vector);
std::ostream& operator<<(std::ostream& os, QuartetHybridRequest const& request);

#endif  // SRC_QUARTET_HYBRID_REQUEST_HPP_
