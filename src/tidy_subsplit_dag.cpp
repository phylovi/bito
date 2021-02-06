// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tidy_subsplit_dag.hpp"

TidySubsplitDAG::TidySubsplitDAG() : SubsplitDAG() {}

TidySubsplitDAG::TidySubsplitDAG(const RootedTreeCollection &tree_collection)
    : TidySubsplitDAG(tree_collection.TaxonCount(), tree_collection.TopologyCounter()) {
}

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

EigenArrayXb TidySubsplitDAG::BelowNode(size_t node_id) {
  return BelowNode(false, node_id).max(BelowNode(true, node_id));
}

EigenArrayXbRef TidySubsplitDAG::BelowNode(bool rotated, size_t node_id) {
  if (rotated) {
    return above_rotated_.col(node_id).array();
  } else {
    return above_sorted_.col(node_id).array();
  }
}

EigenArrayXb TidySubsplitDAG::AboveNode(size_t node_id) {
  return AboveNode(false, node_id).max(AboveNode(true, node_id));
}

EigenArrayXb TidySubsplitDAG::AboveNode(bool rotated, size_t node_id) {
  if (rotated) {
    return above_rotated_.row(node_id).array();
  } else {
    return above_sorted_.row(node_id).array();
  }
}

EigenArrayXbRef TidySubsplitDAG::DirtyVector(bool rotated) {
  if (rotated) {
    return dirty_rotated_;
  } else {
    return dirty_sorted_;
  }
}

bool TidySubsplitDAG::IsDirtyBelow(size_t node_id, bool rotated) {
  // We use `min` as a way of getting "and": we want to find if there are any dirty
  // node-clades below us.
  return BelowNode(rotated, node_id).min(DirtyVector(rotated)).maxCoeff();
}

void TidySubsplitDAG::SetDirtyStrictlyAbove(size_t node_id) {
  for (const bool rotated : {false, true}) {
    EigenArrayXbRef dirty = DirtyVector(rotated);
    EigenArrayXb to_make_dirty = AboveNode(rotated, node_id);
    // We are only dirtying things that are strictly above us.
    to_make_dirty[node_id] = false;
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
  std::stringstream string_stream;
  // I would have thought that we could just do string_stream << m, but this doesn't
  // work.
  for (size_t i = 0; i < m.rows(); i++) {
    string_stream << m.row(i) << "\n";
  }
  return string_stream.str();
}

std::string TidySubsplitDAG::AboveMatricesAsString() {
  std::stringstream string_stream;
  string_stream << "[\n"
                << EigenMatrixXbToString(above_rotated_) << ", \n"
                << EigenMatrixXbToString(above_sorted_) << "\n]";
  return string_stream.str();
}

TidySubsplitDAG TidySubsplitDAG::TrivialExample() {
  // ((0,1),2)
  auto topology = Node::Join(Node::Join(Node::Leaf(0), Node::Leaf(1)), Node::Leaf(2));
  topology->Polish();
  return TidySubsplitDAG(3, {{topology, 1}});
}

TidySubsplitDAG TidySubsplitDAG::ManualTrivialExample() {
  auto manual_dag = TidySubsplitDAG(5);

  // The tree ((0,1)3,2)4:
  // https://github.com/phylovi/libsbn/issues/307#issuecomment-766137769
  manual_dag.SetBelow(3, true, 0);
  manual_dag.SetBelow(3, false, 1);
  manual_dag.SetBelow(4, true, 2);
  manual_dag.SetBelow(4, false, 3);

  return manual_dag;
}

TidySubsplitDAG TidySubsplitDAG::MotivatingExample() {
  auto topologies = Node::ExampleTopologies();
  return TidySubsplitDAG(4, {{topologies[3], 1}, {topologies[4], 1}});
}

std::string TidySubsplitDAG::RecordTraversal() {
  std::stringstream result;
  result << std::boolalpha;
  DepthFirstWithTidyAction(TidySubsplitDAGTraversalAction(
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

void TidySubsplitDAG::SetBelow(size_t parent_id, bool parent_rotated, size_t child_id) {
  BelowNode(parent_rotated, parent_id) =
      BelowNode(parent_rotated, parent_id).max(BelowNode(child_id));
}
