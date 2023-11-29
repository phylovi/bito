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
  // Per-edge adjacent choices for given edge.
  template <typename T>
  struct EdgeAdjacentMap {
    T parent;
    T sister;
    T left_child;
    T right_child;

    T &operator[](EdgeAdjacent edge_adj) {
      switch (edge_adj) {
        case EdgeAdjacent::Parent:
          return parent;
        case EdgeAdjacent::Sister:
          return sister;
        case EdgeAdjacent::LeftChild:
          return left_child;
        case EdgeAdjacent::RightChild:
          return right_child;
      };
      Failwith("ERROR: Invalid enum given.");
    }
    const T &operator[](EdgeAdjacent edge_adj) const {
      switch (edge_adj) {
        case EdgeAdjacent::Parent:
          return parent;
        case EdgeAdjacent::Sister:
          return sister;
        case EdgeAdjacent::LeftChild:
          return left_child;
        case EdgeAdjacent::RightChild:
          return right_child;
      };
      Failwith("ERROR: Invalid enum given.");
    }

    T &operator[](NNIClade clade) {
      switch (clade) {
        case NNIClade::ParentFocal:
          return parent;
        case NNIClade::ParentSister:
          return sister;
        case NNIClade::ChildLeft:
          return left_child;
        case NNIClade::ChildRight:
          return right_child;
      }
      Failwith("ERROR: Invalid enum given.");
    }
    const T &operator[](NNIClade clade) const {
      switch (clade) {
        case NNIClade::ParentFocal:
          return parent;
        case NNIClade::ParentSister:
          return sister;
        case NNIClade::ChildLeft:
          return left_child;
        case NNIClade::ChildRight:
          return right_child;
      }
      Failwith("ERROR: Invalid enum given.");
    }
  };

  using EdgeChoice = EdgeAdjacentMap<EdgeId>;
  using EdgeChoiceNodeIds = EdgeAdjacentMap<NodeId>;
  using EdgeChoiceTreeIds = EdgeAdjacentMap<TreeId>;
  using EdgeChoicePCSPs = NNIAdjacentMap<Bitset>;
  using EdgeChoiceSubsplits = EdgeAdjacentMap<Bitset>;
  using EdgeChoiceNodeIdMap = EdgeAdjacentMap<std::pair<NodeId, NodeId>>;
  using EdgeChoicePCSPMap = EdgeAdjacentMap<std::pair<Bitset, Bitset>>;
  using EdgeChoiceVector = std::vector<EdgeChoice>;
  using TreeIdData = std::vector<TreeId>;

  // ** Constructors

  TPChoiceMap(GPDAG &dag)
      : dag_(dag), edge_choice_vector_(dag.EdgeCountWithLeafSubsplits()){};

  // ** Access

  friend bool operator==(const TPChoiceMap::EdgeChoice &lhs,
                         const TPChoiceMap::EdgeChoice &rhs) {
    if (lhs.parent != rhs.parent) return false;
    if (lhs.sister != rhs.sister) return false;
    if (lhs.left_child != rhs.left_child) return false;
    if (lhs.left_child != rhs.left_child) return false;
    return true;
  }

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
  EdgeId GetEdgeChoice(const EdgeId edge_id, EdgeAdjacent edge_adj) const;
  // Set given edge choice map's given adjacent edge to the given new_edge_choice.
  void SetEdgeChoice(const EdgeId edge_id, const EdgeAdjacent edge_adj,
                     const EdgeId new_edge_choice);
  // Re-initialize edge choices to NoId.
  void ResetEdgeChoice(const EdgeId edge_id);

  // Get data from edge_choice in choice map, according tp edge_id.
  EdgeChoiceNodeIds GetEdgeChoiceNodeIds(const EdgeId edge_id) const;
  EdgeChoicePCSPs GetEdgeChoicePCSPs(const EdgeId edge_id) const;
  EdgeChoiceSubsplits GetEdgeChoiceSubsplits(const EdgeId edge_id) const;

  // Get data from given edge_choice.
  EdgeChoiceNodeIds GetEdgeChoiceNodeIds(const EdgeChoice &edge_choice) const;
  EdgeChoicePCSPs GetEdgeChoicePCSPs(const EdgeChoice &edge_choice) const;
  EdgeChoiceSubsplits GetEdgeChoiceSubsplits(const EdgeChoice &edge_choice) const;

  // Apply NNI clade map transform from pre-NNI to post-NNI.
  template <typename T>
  EdgeAdjacentMap<T> RemapEdgeChoiceDataViaNNICladeMap(
      const EdgeAdjacentMap<T> &old_data,
      const NNIOperation::NNICladeArray &clade_map) const {
    EdgeAdjacentMap<T> new_data{old_data};
    auto SetNewDataFromOldData = [&old_data, &new_data,
                                  &clade_map](const NNIClade nni_clade) {
      new_data[clade_map[nni_clade]] = old_data[nni_clade];
    };
    for (const auto nni_clade : NNICladeEnum::Iterator()) {
      SetNewDataFromOldData(nni_clade);
    }
    return new_data;
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
  // focal edge.
  TreeMask ExtractTreeMask(const EdgeId initial_edge_id) const;
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

  // Extract Topology from DAG with given focal edge.
  Node::Topology ExtractTopology(const EdgeId initial_edge_id) const;
  Node::Topology ExtractTopology(const TreeMask &tree_mask) const;

  // ** I/O

  // Output edge choice map to string by edge_id.
  std::string EdgeChoiceToString(const EdgeId edge_id) const;
  // Output edge choice map to string.
  static std::string EdgeChoiceToString(const EdgeChoice &edge_choice);
  // Output full choice map to string.
  std::string ToString() const;
  // Output choice map as a PCSP Map from the focal edge to vector of adjacent edges.
  using PCSPToPCSPVectorMap = std::map<Bitset, std::vector<Bitset>>;
  PCSPToPCSPVectorMap BuildPCSPMap() const;

  // Output edge choice map.
  friend std::ostream &operator<<(std::ostream &os, const EdgeChoice &edge_choice);
  // Output full choice map.
  friend std::ostream &operator<<(std::ostream &os, const TPChoiceMap &choice_map);

 private:
  // ** ExpandedTreeMask
  // The ExpandedTreeMask contains a map from all the nodes of a TreeMask to their
  // associated parent, left and right child.
  template <typename T>
  using NodeAdjacentArray = EnumArray<NodeAdjacent, 3, T>;
  using ExpandedTreeMask = std::map<NodeId, NodeAdjacentArray<NodeId>>;

  // Extract an ExpandedTreeMask from DAG based on a focal edge or previous TreeMask.
  ExpandedTreeMask ExtractExpandedTreeMask(const EdgeId focal_edge_id) const;
  ExpandedTreeMask ExtractExpandedTreeMask(const TreeMask &tree_mask) const;
  // Extract Tree based on given ExpandedTreeMask.
  Node::Topology ExtractTopology(ExpandedTreeMask &tree_mask_ext) const;
  // Output ExpandedTreeMask to a string.
  std::string ExpandedTreeMaskToString(const ExpandedTreeMask &tree_mask) const;

  // Un-owned reference DAG.
  GPDAG &dag_;
  // A vector that stores a map of each edge's best adjacent edges.
  EdgeChoiceVector edge_choice_vector_;
  // A vector that sets each edge's highest priority tree in the choice map.
  TreeIdData tree_priority_;
};
