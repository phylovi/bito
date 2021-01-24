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
  above_rotated_ = EigenMatrixXb::Identity(node_count, node_count);
  above_sorted_ = EigenMatrixXb::Identity(node_count, node_count);
  dirty_rotated_ = EigenArrayXb::Zero(node_count);
  dirty_sorted_ = EigenArrayXb::Zero(node_count);

  SubsplitDAG::DepthFirstWithAction(SubsplitDAGTraversalAction(
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

EigenArrayXbRef TidySubsplitDAG::DirtyVector(bool rotated) {
  if (rotated) {
    return dirty_rotated_;
  } else {
    return dirty_sorted_;
  }
}

bool TidySubsplitDAG::IsDirtyBelow(size_t node_idx, bool rotated) {
  // We use `min` as a way of getting "and": we want to find if there are any dirty
  // node-clades below us.
  return BelowNode(rotated, node_idx).min(DirtyVector(rotated)).maxCoeff();
}

void TidySubsplitDAG::SetDirtyStrictlyAbove(size_t node_idx) {
  for (const bool rotated : {false, true}) {
    EigenArrayXbRef dirty = DirtyVector(rotated);
    EigenArrayXb to_make_dirty = AboveNode(rotated, node_idx);
    // We are only dirtying things that are strictly above us.
    to_make_dirty[node_idx] = false;
    // We use `max` as a way of getting "or": we want to maintain anything that's
    // already dirty as dirty, while adding all nodes strictly above us.
    dirty = dirty.max(to_make_dirty);
  }
}

void TidySubsplitDAG::SetClean() {
  updating_below_ = std::nullopt;
  dirty_rotated_.setConstant(false);
  dirty_sorted_.setConstant(false);
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

std::string TidySubsplitDAG::RecordTraversal() {
  std::stringstream result;
  result << std::boolalpha;
  DepthFirstWithAction(TidySubsplitDAGTraversalAction(
      // BeforeNode
      [](size_t node_id) {},
      // AfterNode
      [](size_t node_id) {},
      // BeforeNodeClade
      [&result](size_t node_id, bool rotated) {
        result << "descending along " << node_id << ", " << rotated << "\n";
      },
      // ModifyEdge
      [this, &result](size_t node_id, size_t child_id, bool rotated) {
        result << "modifying: ";
        result << node_id << ", " << child_id << ", " << rotated << "\n";
      },
      // UpdateEdge
      [this, &result](size_t node_id, size_t child_id, bool rotated) {
        result << "updating:  ";
        result << node_id << ", " << child_id << ", " << rotated << "\n";
      }));
  return result.str();
}

