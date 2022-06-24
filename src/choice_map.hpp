// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A ChoiceMap is a per-edge map of the best adjacent edges applied to a SubsplitDAG
// for Top-Pruning. Used for selecting, updating, and extracting the top tree from the
// DAG. A ChoiceMap can be used to generate a TreeMask, which is a list of edge ids
// which express a single, complete tree embedded in the DAG, or a Node Topology.

#pragma once

#include <stack>
#include "sugar.hpp"
#include "gp_dag.hpp"
#include "rooted_tree.hpp"

class ChoiceMap {
 public:
  // Per-edge choices of best adjacent edges.
  struct EdgeChoice {
    size_t parent_edge_id = NoId;
    size_t sister_edge_id = NoId;
    size_t left_child_edge_id = NoId;
    size_t right_child_edge_id = NoId;
    double tree_likelihood = -INFINITY;
  };
  using EdgeChoiceVector = std::vector<EdgeChoice>;

  ChoiceMap(GPDAG &dag)
      : dag_(dag), edge_choice_vector_(dag.EdgeCountWithLeafSubsplits()){};

  // ** Selectors

  // Naive choice selector. Chooses the first edge from each list of candidates.
  void SelectFirstEdge();

  // Check if choice selection is valid.
  // Specifically, checks that:
  // - Every edge choice vector has a valid id for all options, unless...
  //   - Edge goes to root (NoId for sister and parent).
  //   - Edge goes to leaf (NoId for left and right child).
  // - Edges span every leaf and root node.
  bool SelectionIsValid(const bool is_quiet = true) const;

  // ** TreeMask

  // A TreeMask is a set of edge Ids, which represent a tree contained in
  // the DAG, from the selected subset of DAG edges.
  using TreeMask = std::set<size_t>;
  // Extract TreeMask from DAG based on edge choices to find best tree with given
  // central edge.
  TreeMask ExtractTreeMask(const size_t central_edge_id) const;
  // Checks that TreeMask represents a valid, complete tree in the DAG.
  // Specifically, checks that:
  // - There is a single edge that goes to the root.
  // - There is a single edge that goes to each leaf.
  // - For each node in mask, there is a single parent, left and right child.
  //   - Unless node is root (no parent) or leaf (no children).
  bool TreeMaskIsValid(const TreeMask &tree_mask, const bool is_quiet = true) const;
  // Output TreeMask to string.
  std::string TreeMaskToString(const TreeMask &tree_mask) const;

  // ** Topology

  // Extract Topology from DAG based on edge choices to find best tree with given
  // central edge.
  Node::NodePtr ExtractTopology(const size_t central_edge_id) const;
  Node::NodePtr ExtractTopology(const TreeMask &tree_mask) const;

  // ** I/O

  // Output edge choice map to iostream.
  friend std::ostream &operator<<(std::ostream &os, const EdgeChoice &edge_choice);
  // Output edge choice map to iostream.
  friend std::ostream &operator<<(std::ostream &os, const ChoiceMap &choice_map);

 private:
  // ** ExpandedTreeMask
  // The ExpandedTreeMask contains a map from all the nodes of a TreeMask to their
  // associated parent, left and right child.
  enum class AdjacentNode { Parent, LeftChild, RightChild };
  template <typename T>
  using AdjacentNodeArray = EnumArray<AdjacentNode, 3, T>;
  using ExpandedTreeMask = std::map<size_t, AdjacentNodeArray<size_t>>;

  // Extract an ExpandedTreeMask from DAG based on a central edge or previous TreeMask.
  ExpandedTreeMask ExtractExpandedTreeMask(const size_t central_edge_id) const;
  ExpandedTreeMask ExtractExpandedTreeMask(const TreeMask &tree_mask) const;
  // Extract Tree based on given ExpandedTreeMask.
  Node::NodePtr ExtractTopology(ExpandedTreeMask &tree_mask_ext) const;
  // Output ExpandedTreeMask to a string.
  std::string ExpandedTreeMaskToString(const ExpandedTreeMask &tree_mask) const;

  // Un-owned reference DAG.
  GPDAG &dag_;
  // A vector that stores a map of each edge's best adjacent edges.
  EdgeChoiceVector edge_choice_vector_;
};
