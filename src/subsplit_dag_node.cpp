// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag_node.hpp"

std::string GetNeighborString(SizeVector neighbors) {
  std::string str;
  for (size_t i : neighbors) {
    str += std::to_string(i) + " ";
  }
  return str;
}

std::string SubsplitDAGNode::ToString() const {
  std::string str = std::to_string(Id()) + ": " + GetBitset().SubsplitToString() + "\n";
  str += "Rootward Sorted: " + GetNeighborString(GetRootwardSorted()) + "\n";
  str += "Rootward Rotated: " + GetNeighborString(GetRootwardRotated()) + "\n";
  str += "Leafward Sorted: " + GetNeighborString(GetLeafwardSorted()) + "\n";
  str += "Leafward Rotated: " + GetNeighborString(GetLeafwardRotated()) + "\n";
  return str;
}
