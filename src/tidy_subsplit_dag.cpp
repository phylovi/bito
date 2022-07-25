// Copyright 2019-2020 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "tidy_subsplit_dag.hpp"

TidySubsplitDAG::TidySubsplitDAG() : SubsplitDAG() {}

TidySubsplitDAG::TidySubsplitDAG(const RootedTreeCollection &tree_collection)
    : TidySubsplitDAG(tree_collection.TaxonCount(), tree_collection.TopologyCounter(),
                      tree_collection.TagTaxonMap()) {}

TidySubsplitDAG::TidySubsplitDAG(size_t node_count)
    : above_rotated_(EigenMatrixXb::Identity(node_count, node_count)),
      above_sorted_(EigenMatrixXb::Identity(node_count, node_count)){};

TidySubsplitDAG::TidySubsplitDAG(size_t taxon_count,
                                 const Node::TopologyCounter &topology_counter,
                                 const TagStringMap &tag_taxon_map)
    : SubsplitDAG(taxon_count, topology_counter, tag_taxon_map) {
  ReinitializeTidyVectors();
}

void TidySubsplitDAG::ReinitializeTidyVectors() {
  auto node_count = NodeCount();
  above_rotated_.resize(node_count, node_count);
  above_rotated_.setIdentity();
  above_sorted_.resize(node_count, node_count);
  above_sorted_.setIdentity();
  dirty_rotated_.resize(node_count);
  dirty_rotated_.setZero();
  dirty_sorted_.resize(node_count);
  dirty_sorted_.setZero();

  SubsplitDAG::DepthFirstWithAction(
      {GetDAGRootNodeId()},
      SubsplitDAGTraversalAction(
          // BeforeNode
          [](NodeId node_id) {},
          // AfterNode
          [](NodeId node_id) {},
          // BeforeNodeClade
          [](NodeId node_id, bool is_edge_on_left) {},
          // VisitEdge
          [this](NodeId node_id, NodeId child_id, bool is_edge_on_left) {
            SetBelow(node_id, is_edge_on_left, child_id);
          }));
}

EigenArrayXb TidySubsplitDAG::BelowNode(NodeId node_id) {
  return BelowNode(false, node_id).max(BelowNode(true, node_id));
}

EigenArrayXbRef TidySubsplitDAG::BelowNode(bool is_edge_on_left, NodeId node_id) {
  if (is_edge_on_left) {
    return above_rotated_.col(node_id.value_).array();
  } else {
    return above_sorted_.col(node_id.value_).array();
  }
}

EigenArrayXb TidySubsplitDAG::AboveNode(NodeId node_id) const {
  return AboveNode(false, node_id).max(AboveNode(true, node_id));
}

EigenArrayXb TidySubsplitDAG::AboveNode(bool is_edge_on_left, NodeId node_id) const {
  if (is_edge_on_left) {
    return above_rotated_.row(node_id.value_).array();
  } else {
    return above_sorted_.row(node_id.value_).array();
  }
}

EigenArrayXbRef TidySubsplitDAG::DirtyVector(bool is_edge_on_left) {
  if (is_edge_on_left) {
    return dirty_rotated_;
  } else {
    return dirty_sorted_;
  }
}

bool TidySubsplitDAG::IsDirtyBelow(NodeId node_id, bool is_edge_on_left) {
  // We use `min` as a way of getting "and": we want to find if there are any dirty
  // node-clades below us.
  return BelowNode(is_edge_on_left, node_id)
      .min(DirtyVector(is_edge_on_left))
      .maxCoeff();
}

void TidySubsplitDAG::SetDirtyStrictlyAbove(NodeId node_id) {
  for (const bool is_edge_on_left : {false, true}) {
    EigenArrayXbRef dirty = DirtyVector(is_edge_on_left);
    EigenArrayXb to_make_dirty = AboveNode(is_edge_on_left, node_id);
    // We are only dirtying things that are strictly above us.
    to_make_dirty[node_id.value_] = false;
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
  for (Eigen::Index i = 0; i < m.rows(); i++) {
    string_stream << m.row(i) << "\n";
  }
  return string_stream.str();
}

std::string TidySubsplitDAG::AboveMatricesAsString() const {
  std::stringstream string_stream;
  string_stream << "[\n"
                << EigenMatrixXbToString(above_rotated_) << ", \n"
                << EigenMatrixXbToString(above_sorted_) << "\n]";
  return string_stream.str();
}

TidySubsplitDAG TidySubsplitDAG::TrivialExample() {
  // ((0,1),2)
  Node::NodePtr topology =
      Node::Join(Node::Join(Node::Leaf(0), Node::Leaf(1)), Node::Leaf(2));
  topology->Polish();
  TagStringMap taxon_map = TidySubsplitDAG::BuildDummyTagTaxonMap(3);
  return TidySubsplitDAG(3, {{topology, 1}}, taxon_map);
}

TidySubsplitDAG TidySubsplitDAG::ManualTrivialExample() {
  auto manual_dag = TidySubsplitDAG(6);

  // The tree ((0,1)3,2)4:
  // https://github.com/phylovi/bito/issues/349#issuecomment-897963382
  manual_dag.SetBelow(NodeId(3), true, NodeId(0));
  manual_dag.SetBelow(NodeId(3), false, NodeId(1));
  manual_dag.SetBelow(NodeId(4), false, NodeId(2));
  manual_dag.SetBelow(NodeId(4), true, NodeId(3));
  manual_dag.SetBelow(NodeId(5), true, NodeId(4));

  return manual_dag;
}

TidySubsplitDAG TidySubsplitDAG::MotivatingExample() {
  auto topologies = Node::ExampleTopologies();
  TagStringMap taxon_map = TidySubsplitDAG::BuildDummyTagTaxonMap(4);
  return TidySubsplitDAG(4, {{topologies[3], 1}, {topologies[4], 1}}, taxon_map);
}

std::string TidySubsplitDAG::RecordTraversal() {
  std::stringstream result;
  result << std::boolalpha;
  DepthFirstWithTidyAction(
      {GetDAGRootNodeId()},
      TidySubsplitDAGTraversalAction(
          // BeforeNode
          [](NodeId node_id) {},
          // AfterNode
          [](NodeId node_id) {},
          // BeforeNodeClade
          [&result](NodeId node_id, bool is_edge_on_left) {
            result << "descending along " << node_id << ", " << is_edge_on_left << "\n";
          },
          // ModifyEdge
          [&result](NodeId node_id, NodeId child_id, bool is_edge_on_left) {
            result << "modifying: ";
            result << node_id << ", " << child_id << ", " << is_edge_on_left << "\n";
          },
          // UpdateEdge
          [&result](NodeId node_id, NodeId child_id, bool is_edge_on_left) {
            result << "updating:  ";
            result << node_id << ", " << child_id << ", " << is_edge_on_left << "\n";
          }));
  return result.str();
}

void TidySubsplitDAG::SetBelow(NodeId parent_id, bool is_edge_on_left,
                               NodeId child_id) {
  BelowNode(is_edge_on_left, parent_id) =
      BelowNode(is_edge_on_left, parent_id).max(BelowNode(child_id));
}
