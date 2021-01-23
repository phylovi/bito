// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tidy_subsplit_dag.hpp"

TidySubsplitDAG::TidySubsplitDAG() : SubsplitDAG() {}

// TODO something
TidySubsplitDAG::TidySubsplitDAG(const RootedTreeCollection &tree_collection)
    : SubsplitDAG(tree_collection) {}

TidySubsplitDAG::TidySubsplitDAG(size_t node_count)
    : above_sorted_(EigenMatrixXb::Identity(node_count, node_count)),
      above_rotated_(EigenMatrixXb::Identity(node_count, node_count)){};

TidySubsplitDAG::TidySubsplitDAG(size_t taxon_count,
                                 const Node::TopologyCounter &topology_counter)
    : SubsplitDAG(taxon_count, topology_counter) {
  auto node_count = NodeCount();
  above_sorted_ = EigenMatrixXb::Identity(node_count, node_count);
  above_rotated_ = EigenMatrixXb::Identity(node_count, node_count);

  DepthFirstWithAction(SubsplitDAGTraversalAction(
      // BeforeNode
      [](size_t node_id) {},
      // AfterNode
      [](size_t node_id) {},
      // BeforeNodeClade
      [](size_t node_id, bool rotated) {},
      // VisitEdge
      [this](size_t node_id, size_t child_id, bool rotated) {
        SetBelow(node_id, rotated, child_id);
      }));
}

EigenArrayXb TidySubsplitDAG::BelowNode(size_t node_idx) {
  return BelowNode(false, node_idx).max(BelowNode(true, node_idx));
}

EigenArrayXbRef TidySubsplitDAG::BelowNode(bool rotated, size_t node_idx) {
  if (rotated) {
    return above_rotated_.col(node_idx).array();
  } else {
    return above_sorted_.col(node_idx).array();
  }
}

EigenArrayXb TidySubsplitDAG::AboveNode(size_t node_idx) {
  return AboveNode(false, node_idx).max(AboveNode(true, node_idx));
}

EigenArrayXb TidySubsplitDAG::AboveNode(bool rotated, size_t node_idx) {
  if (rotated) {
    return above_rotated_.row(node_idx).array();
  } else {
    return above_sorted_.row(node_idx).array();
  }
}

// TODO change idx to id
void TidySubsplitDAG::SetBelow(size_t parent_idx, bool parent_rotated,
                               size_t child_idx) {
  BelowNode(parent_rotated, parent_idx) =
      BelowNode(parent_rotated, parent_idx).max(BelowNode(child_idx));
}

std::string EigenMatrixXbToString(EigenMatrixXb m) {
  std::stringstream ss;
  // I would have thought that we could just do ss << m, but this doesn't work.
  for (size_t i = 0; i < m.rows(); i++) {
    ss << m.row(i) << "\n";
  }
  return ss.str();
}

std::string TidySubsplitDAG::AboveMatricesAsString() {
  std::stringstream ss;
  ss << "[\n"
     << EigenMatrixXbToString(above_rotated_) << ", \n"
     << EigenMatrixXbToString(above_sorted_) << "\n]";
  return ss.str();
}

TidySubsplitDAG TidySubsplitDAG::TrivialExample() {
  // ((0,1),2)
  auto topology = Node::Join(Node::Join(Node::Leaf(0), Node::Leaf(1)), Node::Leaf(2));
  topology->Polish();
  return TidySubsplitDAG(3, {{topology, 1}});
}

TidySubsplitDAG TidySubsplitDAG::MotivatingExample() {
  auto topologies = Node::ExampleTopologies();
  return TidySubsplitDAG(4, {{topologies[3], 1}, {topologies[4], 1}});
}

