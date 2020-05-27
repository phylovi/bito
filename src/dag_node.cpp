// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "dag_node.hpp"

std::string GetNeighborString(std::vector<size_t> neighbors)
{
  std::string str = "";
  for (size_t i : neighbors) {
    str += std::to_string(i) + " ";
  }
  return str;
}

bool DAGNode::IsLeaf() const {
  return (leafward_rotated.size() == 0) && (leafward_sorted.size() == 0);
}

void DAGNode::AddNeighbor(EdgeType edge_type, size_t node_id)
{
  switch (edge_type) {
    case EdgeType::LEAFWARD_ROTATED:
      leafward_rotated.push_back(node_id);
      break;
    case EdgeType::LEAFWARD_SORTED:
      leafward_sorted.push_back(node_id);
      break;
    case EdgeType::ROOTWARD_SORTED:
      rootward_sorted.push_back(node_id);
      break;
    case EdgeType::ROOTWARD_ROTATED:
      rootward_rotated.push_back(node_id);
      break;
  }
}

std::string DAGNode::ToString()
{
  std::string str = std::to_string(id_) + "\n";
  str += "Rootward Sorted: " + GetNeighborString(rootward_sorted) + "\n";
  str += "Rootward Rotated: " + GetNeighborString(rootward_rotated) + "\n";
  str += "Leafward Sorted: " + GetNeighborString(leafward_sorted) + "\n";
  str += "Leafward Rotated: " + GetNeighborString(leafward_rotated) + "\n";
  return str;
}
