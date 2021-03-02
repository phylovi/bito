// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "quartet_hybrid_request.hpp"

bool QuartetHybridRequest::IsComplete() const {
  return !rootward_tips_.empty() && !sister_tips_.empty() && !rotated_tips_.empty() &&
         !sorted_tips_.empty();
}

std::ostream& operator<<(std::ostream& os, QuartetTip const& plv_pcsp) {
  os << "(tip node " << plv_pcsp.tip_node_id_ << ", plv " << plv_pcsp.plv_idx_
     << ", gpcsp " << plv_pcsp.gpcsp_idx_ << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, QuartetTipVector const& plv_pcsp_vector) {
  os << "[";
  for (const auto& plv_pcsp : plv_pcsp_vector) {
    os << plv_pcsp << ", ";
  }
  os << "]";
  return os;
}

std::ostream& operator<<(std::ostream& os, QuartetHybridRequest const& request) {
  os << "[\n";
  os << "\tcentral GPCSP: " << request.central_gpcsp_idx_ << "\n";
  os << "\trootward tips: " << request.rootward_tips_ << "\n";
  os << "\tsister tips: " << request.sister_tips_ << "\n";
  os << "\trotated tips: " << request.rotated_tips_ << "\n";
  os << "\tsorted tips: " << request.sorted_tips_ << "\n";
  os << "]" << std::endl;
  return os;
}
