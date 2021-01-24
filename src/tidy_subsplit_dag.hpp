// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A "tidy" subsplit DAG has a notion of clean and dirty vectors.

#ifndef SRC_TIDY_SUBSPLIT_DAG_HPP_
#define SRC_TIDY_SUBSPLIT_DAG_HPP_

#include "subsplit_dag.hpp"

class TidySubsplitDAG : public SubsplitDAG {
 public:
  TidySubsplitDAG();
  explicit TidySubsplitDAG(const RootedTreeCollection &tree_collection);
  // This constructor is really just meant for testing.
  explicit TidySubsplitDAG(size_t node_count);

  // What nodes are above or below the specified node? We consider a node to be both
  // above and below itself (this just happens to be handy for the implementation).
  EigenArrayXb BelowNode(size_t node_idx);
  // TODO let's change to rotated second. It's a subsplit clade.
  EigenArrayXbRef BelowNode(bool rotated, size_t node_idx);
  EigenArrayXb AboveNode(size_t node_idx);
  EigenArrayXb AboveNode(bool rotated, size_t node_idx);

  // Set the below matrix up to have the node at src_idx below the subsplit-clade
  // described by (dst_rotated, dst_idx).
  void SetBelow(size_t dst_idx, bool dst_rotated, size_t src_idx);

  EigenArrayXbRef DirtyVector(bool rotated);
  bool IsDirtyBelow(size_t node_idx, bool rotated);
  void SetDirtyAbove(size_t node_idx, bool rotated);

  std::string AboveMatricesAsString();

  // From ((0,1)2)
  // https://github.com/phylovi/libsbn/issues/307#issuecomment-766137769
  static TidySubsplitDAG TrivialExample();
  // From (0,(1,(2,3))) and ((0,(2,3)),1)
  // See https://github.com/phylovi/libsbn/issues/307#issuecomment-765901588
  // Update during #288
  static TidySubsplitDAG MotivatingExample();

  /*
  // Apply a TidySubsplitDAGTraversalAction via a depth first traversal. Do not visit
  // leaf nodes. Applied to a given node, we:
  // * Apply BeforeNode
  // * For each of the clades of the node, we:
  //     * Descend into each clade, cleaning up with UpdateEdge as needed.
  //     * Apply BeforeNodeClade
  //     * For each edge descending from that clade, we:
  //         * Recur into the child node of the clade if it is not a leaf
  //         * Apply VisitEdge to the edge
  // * Apply AfterNode
  template <typename TraversalActionT>
  void DepthFirstWithAction(const TraversalActionT &action) const {
    std::unordered_set<size_t> visited_nodes;
    for (const auto &rootsplit : rootsplits_) {
      DepthFirstWithActionForNode(action, subsplit_to_id_.at(rootsplit + ~rootsplit),
                                  visited_nodes);
    }
  };

  // The portion of the traversal that is below a given node.
  template <typename TraversalActionT>
  void DepthFirstWithActionForNode(const TraversalActionT &action, size_t node_id,
                                   std::unordered_set<size_t> &visited_nodes) const {
    action.BeforeNode(node_id);
    DepthFirstWithActionForNodeClade(action, node_id, false, visited_nodes);
    DepthFirstWithActionForNodeClade(action, node_id, true, visited_nodes);
    action.AfterNode(node_id);
  };

  // The portion of the traversal that is below a given clade of a given node.
  // Do not recur into leaf nodes.
  template <typename TraversalActionT>
  void DepthFirstWithActionForNodeClade(
      const TraversalActionT &action, size_t node_id, bool rotated,
      std::unordered_set<size_t> &visited_nodes) const {
    // TODO clean up with UpdateEdge
    action.BeforeNodeClade(node_id, rotated);
    const auto node = GetDagNode(node_id);
    for (const size_t child_id : node->GetLeafward(rotated)) {
      if (visited_nodes.count(child_id) == 0) {
        visited_nodes.insert(child_id);
        if (!GetDagNode(child_id)->IsLeaf()) {
          DepthFirstWithActionForNode(action, child_id, visited_nodes);
        }
      }
      action.ModifyEdge(node_id, child_id, rotated);
      // TODO dirty edges above us
    }
  };
  */

 private:
  // If this is set then we are in an "updating mode", where we are updating below the
  // specified node-clade.
  std::optional<std::pair<bool, size_t>> updating_below;

  // above_rotated_(i,j) is true iff i,true is above j.
  EigenMatrixXb above_rotated_;
  // above_sorted(i,j) is true iff i,false is above j.
  EigenMatrixXb above_sorted_;

  // dirty_rotated_(i) is true iff i,true is dirty.
  EigenArrayXb dirty_rotated_;
  // dirty_sorted_(i) is true iff i,false is dirty.
  EigenArrayXb dirty_sorted_;

  TidySubsplitDAG(size_t taxon_count, const Node::TopologyCounter &topology_counter);
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TidySubsplitDAG: slicing") {
  auto manual_dag = TidySubsplitDAG(5);

  // The tree ((0,1)3,2)4:
  // https://github.com/phylovi/libsbn/issues/307#issuecomment-766137769
  manual_dag.SetBelow(3, true, 0);
  manual_dag.SetBelow(3, false, 1);
  manual_dag.SetBelow(4, true, 2);
  manual_dag.SetBelow(4, false, 3);

  CHECK_EQ(GenericToString(manual_dag.AboveNode(0)), "[1, 0, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(manual_dag.AboveNode(1)), "[0, 1, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(manual_dag.AboveNode(2)), "[0, 0, 1, 0, 1]\n");
  CHECK_EQ(GenericToString(manual_dag.AboveNode(3)), "[0, 0, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(manual_dag.AboveNode(4)), "[0, 0, 0, 0, 1]\n");

  auto trivial_dag = TidySubsplitDAG::TrivialExample();
  CHECK_EQ(trivial_dag.AboveMatricesAsString(), manual_dag.AboveMatricesAsString());

  auto motivating_dag = TidySubsplitDAG::MotivatingExample();
  CHECK_EQ(GenericToString(motivating_dag.AboveNode(false, 4)),
           "[0, 0, 0, 0, 1, 1, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(motivating_dag.AboveNode(true, 4)),
           "[0, 0, 0, 0, 1, 0, 1, 0, 0]\n");
  CHECK_EQ(GenericToString(motivating_dag.AboveNode(false, 7)),
           "[0, 0, 0, 0, 0, 0, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(motivating_dag.AboveNode(true, 7)),
           "[0, 0, 0, 0, 0, 0, 0, 1, 0]\n");
  CHECK_EQ(GenericToString(motivating_dag.BelowNode(false, 7)),
           "[0, 0, 1, 1, 1, 0, 0, 1, 0]\n");
  CHECK_EQ(GenericToString(motivating_dag.BelowNode(true, 7)),
           "[1, 0, 0, 0, 0, 0, 0, 1, 0]\n");
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TIDY_SUBSPLIT_DAG_HPP_
