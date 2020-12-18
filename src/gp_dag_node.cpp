// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_dag_node.hpp"

std::string GetNeighborString(std::vector<size_t> neighbors) {
  std::string str;
  for (size_t i : neighbors) {
    str += std::to_string(i) + " ";
  }
  return str;
}

std::string GPDAGNode::ToString() {
  std::string str = std::to_string(id_) + ": " + GetBitset().SubsplitToString() + "\n";
  str += "Rootward Sorted: " + GetNeighborString(rootward_sorted_) + "\n";
  str += "Rootward Rotated: " + GetNeighborString(rootward_rotated_) + "\n";
  str += "Leafward Sorted: " + GetNeighborString(leafward_sorted_) + "\n";
  str += "Leafward Rotated: " + GetNeighborString(leafward_rotated_) + "\n";
  return str;
}
