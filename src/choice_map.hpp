// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A ChoiceMap is a per-edge map of the best adjacent edges applied to a SubsplitDAG
// for Top-Pruning. Used for selecting, updating, and extracting the top tree from the
// DAG. A ChoiceMap can generate a TreeMask, which is a list of edge ids which express a
// single, complete tree embedded in the DAG.

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
  enum class AdjacentNode { Parent, LeftChild, RightChild };
  template <typename T>
  using AdjacentNodeArray = EnumArray<AdjacentNode, 3, T>;
  using ExpandedTreeMask = std::map<size_t, AdjacentNodeArray<size_t>>;

  // Extract TreeMask from DAG based on edge choices to find best tree with given
  // central edge.
  // - Makes two passes:
  //   - The first pass goes up along the chosen edges of the DAG to the root, adding
  //   each edge it encounters.
  //   - The second pass goes leafward, descending along the chosen edges to the leaf
  //   edges from the sister of each edge in the rootward pass and the child edges from
  //   the central edge.
  TreeMask ExtractTreeMask(const size_t central_edge_id) const;

  // Expand a TreeMask into a ExpandedTreeMask.  Rather than containing all the edges in
  // the tree, the ExpandedTreeMask contains all the nodes in the tree in a map
  // associated with their parent, left and right child.
  ExpandedTreeMask ExtractExpandedTreeMask(const size_t central_edge_id) const;
  ExpandedTreeMask ExtractExpandedTreeMask(const TreeMask &tree_mask) const;

  // Checks that TreeMask represents a valid, complete tree in the DAG.
  // Specifically, checks that:
  // - There is a single edge that goes to the root.
  // - There is a single edge that goes to each leaf.
  // - For each node in mask, there is a single parent, left and right child.
  //   - Unless node is root (no parent) or leaf (no children).
  bool TreeMaskIsValid(const TreeMask &tree_mask, const bool is_quiet = true) const;

  // Output TreeMask to string.
  std::string TreeMaskToString(const TreeMask &tree_mask) const;
  std::string ExpandedTreeMaskToString(const ExpandedTreeMask &tree_mask) const;

  // ** Tree

  // Extract Tree from DAG based on edge choices to find best tree with given
  // central edge.
  // - Makes two passes:
  //   - The first pass goes up along the chosen edges of the DAG to the root, adding
  //   each edge it encounters.
  //   - The second pass goes leafward, descending along the chosen edges to the leaf
  //   edges from the sister of each edge in the rootward pass and the child edges from
  //   the central edge.
  Tree ExtractTree(const size_t central_edge_id) const;
  Tree ExtractTree(const TreeMask &tree_mask) const;
  Tree ExtractTree(ExpandedTreeMask &tree_mask_ext) const;

  // Checks that tree is a valid tree in DAG that spans the
  bool TreeIsValid(const Tree &tree, const bool is_quiet = true) const;

  // ** I/O

  // Output edge choice map to iostream.
  friend std::ostream &operator<<(std::ostream &os, const EdgeChoice &edge_choice);
  // Output edge choice map to iostream.
  friend std::ostream &operator<<(std::ostream &os, const ChoiceMap &choice_map);

 private:
  // Un-owned reference DAG.
  GPDAG &dag_;
  // A vector that stores a map of each edge's best adjacent edges.
  EdgeChoiceVector edge_choice_vector_;
};
