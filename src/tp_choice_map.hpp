// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A TPChoiceMap is a per-edge map of the best adjacent edges applied to a SubsplitDAG
// for Top-Pruning. Used for selecting, updating, and extracting the top tree from the
// DAG. A TPChoiceMap can be used to generate a TreeMask, which is a list of edge ids
// which express a single, complete tree embedded in the DAG, or a Node Topology.

#pragma once

#include <stack>
#include "sugar.hpp"
#include "gp_dag.hpp"
#include "node.hpp"
#include "dag_data.hpp"

using TreeId = GenericId<struct TreeIdTag>;

class TPChoiceMap {
 public:
  // Types of Adjacent Nodes
  enum class AdjacentNode { Parent, LeftChild, RightChild };
  static const inline size_t AdjacentNodeCount = 3;
  class AdjacentNodeEnum
      : public EnumWrapper<AdjacentNode, size_t, AdjacentNodeCount,
                           AdjacentNode::Parent, AdjacentNode::RightChild> {
   public:
    static inline const std::string Prefix = "AdjacentNode";
    static inline const Array<std::string> Labels = {
        {"Parent", "LeftChild", "RightChild"}};

    static std::string ToString(const AdjacentNode e) {
      std::stringstream ss;
      ss << Prefix << "::" << Labels[e];
      return ss.str();
    }
    friend std::ostream &operator<<(std::ostream &os, const AdjacentNode e) {
      os << ToString(e);
      return os;
    }
  };

  // Types of Adjacent Edges
  enum class AdjacentEdge { Parent, Sister, LeftChild, RightChild };
  static const inline size_t AdjacentEdgeCount = 4;
  class AdjacentEdgeEnum
      : public EnumWrapper<AdjacentEdge, size_t, AdjacentEdgeCount,
                           AdjacentEdge::Parent, AdjacentEdge::RightChild> {
   public:
    static inline const std::string Prefix = "AdjacentEdge";
    static inline const Array<std::string> Labels = {
        {"Parent", "Sister", "LeftChild", "RightChild"}};

    static std::string ToString(const AdjacentEdge e) {
      std::stringstream ss;
      ss << Prefix << "::" << Labels[e];
      return ss.str();
    }
    friend std::ostream &operator<<(std::ostream &os, const AdjacentEdge e) {
      os << ToString(e);
      return os;
    }
  };

  // Per-edge choices of best adjacent edges.
  struct EdgeChoice {
    EdgeId parent_edge_id = EdgeId(NoId);
    EdgeId sister_edge_id = EdgeId(NoId);
    EdgeId left_child_edge_id = EdgeId(NoId);
    EdgeId right_child_edge_id = EdgeId(NoId);
  };
  using EdgeChoiceVector = std::vector<EdgeChoice>;

  friend bool operator==(const TPChoiceMap::EdgeChoice &lhs,
                         const TPChoiceMap::EdgeChoice &rhs) {
    if (lhs.parent_edge_id != rhs.parent_edge_id) return false;
    if (lhs.sister_edge_id != rhs.sister_edge_id) return false;
    if (lhs.left_child_edge_id != rhs.left_child_edge_id) return false;
    if (lhs.left_child_edge_id != rhs.left_child_edge_id) return false;
    return true;
  }

  // Per-edge choices of best adjacent nodes. Note: these are the nearest common nodes
  // of an NNI Operation.
  struct EdgeChoiceNodeIds {
    NodeId parent_node_id = NodeId(NoId);
    NodeId sister_node_id = NodeId(NoId);
    NodeId left_child_node_id = NodeId(NoId);
    NodeId right_child_node_id = NodeId(NoId);
  };

  // Per-edge choices associated tree_ids.
  struct EdgeChoiceTreeIds {
    TreeId parent_tree_id = TreeId(NoId);
    TreeId sister_tree_id = TreeId(NoId);
    TreeId left_child_tree_id = TreeId(NoId);
    TreeId right_child_tree_id = TreeId(NoId);
  };

  struct EdgeChoicePCSPs {
    Bitset parent_pcsp = Bitset(0);
    Bitset focal_pcsp = Bitset(0);
    Bitset sister_pcsp = Bitset(0);
    Bitset left_child_pcsp = Bitset(0);
    Bitset right_child_pcsp = Bitset(0);
  };

  TPChoiceMap(GPDAG &dag)
      : dag_(dag), edge_choice_vector_(dag.EdgeCountWithLeafSubsplits()){};

  // ** Access

  // Size of edge choice map.
  size_t size() const { return edge_choice_vector_.size(); }
  // Get associated DAG.
  const GPDAG &GetDAG() const { return dag_; }

  // Get choice map for given edge_id.
  EdgeChoice &GetEdgeChoice(const EdgeId edge_id) {
    return edge_choice_vector_[edge_id.value_];
  }
  const EdgeChoice &GetEdgeChoice(const EdgeId edge_id) const {
    return edge_choice_vector_[edge_id.value_];
  }

  // Get adjacent edge id in given edge's choice map for adjacent edge direction.
  EdgeId GetEdgeChoice(const EdgeId edge_id, AdjacentEdge edge_choice_type) {
    switch (edge_choice_type) {
      case AdjacentEdge::Parent:
        return edge_choice_vector_[edge_id.value_].parent_edge_id;
      case AdjacentEdge::Sister:
        return edge_choice_vector_[edge_id.value_].sister_edge_id;
      case AdjacentEdge::LeftChild:
        return edge_choice_vector_[edge_id.value_].left_child_edge_id;
      case AdjacentEdge::RightChild:
        return edge_choice_vector_[edge_id.value_].right_child_edge_id;
      default:
        Failwith("Invalid edge choice type.");
    }
  }

  // Set given edge choice map's given adjacent edge to the given new_edge_choice.
  void SetEdgeChoice(const EdgeId edge_id, const AdjacentEdge edge_choice_type,
                     const EdgeId new_edge_choice) {
    switch (edge_choice_type) {
      case AdjacentEdge::Parent:
        edge_choice_vector_[edge_id.value_].parent_edge_id = new_edge_choice;
        break;
      case AdjacentEdge::Sister:
        edge_choice_vector_[edge_id.value_].sister_edge_id = new_edge_choice;
        break;
      case AdjacentEdge::LeftChild:
        edge_choice_vector_[edge_id.value_].left_child_edge_id = new_edge_choice;
        break;
      case AdjacentEdge::RightChild:
        edge_choice_vector_[edge_id.value_].right_child_edge_id = new_edge_choice;
        break;
      default:
        Failwith("Invalid edge choice type.");
    }
  }
  // Re-initialize edge choices to NoId.
  void ResetEdgeChoice(const EdgeId edge_id) {
    edge_choice_vector_[edge_id.value_] = {EdgeId(NoId), EdgeId(NoId), EdgeId(NoId),
                                           EdgeId(NoId)};
  }

  // Get the node ids corresponding to a given edge_id's choice map.
  EdgeChoiceNodeIds GetEdgeChoiceNodeIds(const EdgeId edge_id) const {
    const auto &edge_choice = GetEdgeChoice(edge_id);
    return GetNodeIdsFromEdgeChoice(edge_choice);
  }
  EdgeChoiceNodeIds GetNodeIdsFromEdgeChoice(const EdgeChoice &edge_choice) const {
    auto GetNodeFromEdge = [this](const EdgeId edge_id, const Direction direction) {
      if (edge_id == NoId) {
        return NodeId(NoId);
      }
      const auto &edge = GetDAG().GetDAGEdge(edge_id);
      return (direction == Direction::Rootward) ? edge.GetParent() : edge.GetChild();
    };

    EdgeChoiceNodeIds node_choice;
    node_choice.parent_node_id =
        GetNodeFromEdge(edge_choice.parent_edge_id, Direction::Rootward);
    node_choice.sister_node_id =
        GetNodeFromEdge(edge_choice.sister_edge_id, Direction::Leafward);
    node_choice.left_child_node_id =
        GetNodeFromEdge(edge_choice.left_child_edge_id, Direction::Leafward);
    node_choice.right_child_node_id =
        GetNodeFromEdge(edge_choice.right_child_edge_id, Direction::Leafward);
    return node_choice;
  }

  EdgeChoicePCSPs GetEdgeChoicePCSPs(const EdgeId edge_id) const {
    EdgeChoicePCSPs adj_pcsps;
    const auto &edge_choice = GetEdgeChoice(edge_id);
    auto SetPCSPIfExists = [this](Bitset &adj_edge_bitset, const EdgeId adj_edge_id) {
      if (adj_edge_id != NoId) {
        adj_edge_bitset = GetDAG().GetDAGEdgeBitset(adj_edge_id);
      }
    };
    SetPCSPIfExists(adj_pcsps.parent_pcsp, edge_choice.parent_edge_id);
    SetPCSPIfExists(adj_pcsps.focal_pcsp, edge_id);
    SetPCSPIfExists(adj_pcsps.sister_pcsp, edge_choice.sister_edge_id);
    SetPCSPIfExists(adj_pcsps.left_child_pcsp, edge_choice.left_child_edge_id);
    SetPCSPIfExists(adj_pcsps.right_child_pcsp, edge_choice.right_child_edge_id);

    return adj_pcsps;
  }

  // ** Maintenance

  // Grow and reindex data to fit new DAG. Initialize new choice map to first edge.
  void GrowEdgeData(const size_t new_edge_count,
                    std::optional<const Reindexer> edge_reindexer,
                    std::optional<const size_t> explicit_alloc, const bool on_init);

  // ** Selectors

  // Naive choice selector. Chooses the first edge from each list of candidates.
  void SelectFirstEdge();
  void SelectFirstEdge(const EdgeId edge_id);

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
  using TreeMask = std::set<EdgeId>;
  // Extract TreeMask from DAG based on edge choices to find best tree with given
  // central edge.
  TreeMask ExtractTreeMask(const EdgeId central_edge_id) const;
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
  // Extract tree topology from DAG based on edges choices to find best tree.

  // Extract Topology from DAG with given central edge.
  Node::Topology ExtractTopology(const EdgeId central_edge_id) const;
  Node::Topology ExtractTopology(const TreeMask &tree_mask) const;

  // ** I/O

  // Output edge choice map to string by edge_id.
  std::string EdgeChoiceToString(const EdgeId edge_id) const;
  // Output edge choice map to string.
  static std::string EdgeChoiceToString(const EdgeChoice &edge_choice);
  // Output full choice map to string.
  std::string ToString() const;

  // Output edge choice map.
  friend std::ostream &operator<<(std::ostream &os, const EdgeChoice &edge_choice);
  // Output full choice map.
  friend std::ostream &operator<<(std::ostream &os, const TPChoiceMap &choice_map);

 private:
  // ** ExpandedTreeMask
  // The ExpandedTreeMask contains a map from all the nodes of a TreeMask to their
  // associated parent, left and right child.
  template <typename T>
  using AdjacentNodeArray = EnumArray<AdjacentNode, 3, T>;
  using ExpandedTreeMask = std::map<NodeId, AdjacentNodeArray<NodeId>>;

  // Extract an ExpandedTreeMask from DAG based on a central edge or previous TreeMask.
  ExpandedTreeMask ExtractExpandedTreeMask(const EdgeId central_edge_id) const;
  ExpandedTreeMask ExtractExpandedTreeMask(const TreeMask &tree_mask) const;
  // Extract Tree based on given ExpandedTreeMask.
  Node::Topology ExtractTopology(ExpandedTreeMask &tree_mask_ext) const;
  // Output ExpandedTreeMask to a string.
  std::string ExpandedTreeMaskToString(const ExpandedTreeMask &tree_mask) const;

  // Un-owned reference DAG.
  GPDAG &dag_;
  // A vector that stores a map of each edge's best adjacent edges.
  EdgeChoiceVector edge_choice_vector_;
  // A vector that sets the priority of each tree in the choice map.
  DAGEdgeIntData tree_priority_;
};
