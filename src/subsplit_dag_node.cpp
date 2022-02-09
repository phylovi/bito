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

bool SubsplitDAGNode::IsValid() const {
  // If node is a leaf, then a valid node should have no parents.
  if (IsLeaf()) {
    return (GetRightLeafward().size() + GetLeftLeafward().size() == 0);
  }
  // If node is a root, then a valid node should have no children.
  else if (IsDAGRootNode()) {
    return (GetRightRootward().size() + GetLeftRootward().size() == 0);
  }
  // If neither, then node should either have:
  // (1) Zero parents and zero children.
  // (2) 1+ parents, 1+ sorted children, and 1+ rotated children.
  size_t parent_node_count = GetRightRootward().size() + GetLeftRootward().size();
  if (parent_node_count > 0) {
    if (GetRightLeafward().size() == 0 || GetRightLeafward().size() == 0) {
      return false;
    }
  } else {
    if (GetRightLeafward().size() > 0 || GetRightLeafward().size() > 0) {
      return false;
    }
  }
  return true;
}

std::string SubsplitDAGNode::ToString() const {
  std::string str = std::to_string(id_) + ": " + GetBitset().SubsplitToString() + "\n";
  str += "Right Rootward: " + GetNeighborString(right_rootward_) + "\n";
  str += "Left Rootward: " + GetNeighborString(left_rootward_) + "\n";
  str += "Right Leafward: " + GetNeighborString(right_leafward_) + "\n";
  str += "Left Leafward: " + GetNeighborString(left_leafward_) + "\n";
  return str;
}
