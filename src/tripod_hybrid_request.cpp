// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tripod_hybrid_request.hpp"

std::ostream& operator<<(std::ostream& os, TripodTip const& plv_pcsp) {
  os << "(tip node " << plv_pcsp.tip_node_id_ << ", plv " << plv_pcsp.plv_idx_
     << ", gpcsp " << plv_pcsp.gpcsp_idx_ << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, TripodTipVector const& plv_pcsp_vector) {
  os << "[";
  for (const auto& plv_pcsp : plv_pcsp_vector) {
    os << plv_pcsp << ", ";
  }
  os << "]";
  return os;
}

std::ostream& operator<<(std::ostream& os, TripodHybridRequest const& request) {
  os << "[\n";
  os << "\tcentral GPCSP: " << request.central_gpcsp_idx_ << "\n";
  os << "\trootward pairs: " << request.rootward_tips_ << "\n";
  os << "\trotated pairs: " << request.rotated_tips_ << "\n";
  os << "\tsorted pairs: " << request.sorted_tips_ << "\n";
  os << "]" << std::endl;
  return os;
}
