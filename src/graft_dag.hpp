// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The GraftDAG is a proposed addition (graft) to SubsplitDAG (host), which we
// can perform computations on without the need for adding nodes and edges and
// reindexing the full DAG.

#include "gp_dag.hpp"
#include "subsplit_dag.hpp"

#pragma once

class GraftDAG {
 public:
  // ** Constructors:

  // Initialize empty GraftDAG.
  GraftDAG(SubsplitDAG &dag);
  // Initialize GraftDAG with an initial graft.
  GraftDAG(SubsplitDAG &dag, NNIOperation &nni);
  GraftDAG(SubsplitDAG &dag, Bitset &parent_subsplit, Bitset &child_subsplit);

  // ** Comparators

  // Uses same method of comparison as SubsplitDAG (node and edge sets).
  int Compare(const GraftDAG &other) const;
  static int Compare(const GraftDAG &lhs, const GraftDAG &rhs);
  // Treats GraftDAG as completed DAG to compare against normal SubsplitDAG.
  int CompareToDAG(const SubsplitDAG &other) const;
  static int CompareToDAG(const GraftDAG &lhs, const SubsplitDAG &rhs);

  // ** Modify GraftDAG

  // Add node pair to graft.
  void AddGraftNodePair(const NNIOperation &nni);
  void AddGraftNodePair(const Bitset &parent_subsplit, const Bitset &child_subsplit);

  // Clear all nodes and edges from graft for reuse.
  void RemoveAllGrafts();

  // ** Build Indexers/Vectors

  // Get the left and right parents of the node with the given subsplit.
  std::pair<SizeVector, SizeVector> BuildParentIdVectors(const Bitset &subsplit) const;
  // Get the left and right children of the node with the given subsplit.
  std::pair<SizeVector, SizeVector> BuildChildIdVectors(const Bitset &subsplit) const;

  // ** Getters

  // Get pointer to the host DAG.
  const SubsplitDAG &GetHostDAG() const;
  // Get node based on node id.
  SubsplitDAGNode GetDAGNode(const size_t node_id) const;
  MutableSubsplitDAGNode GetDAGNode(const size_t node_id);
  // Get the node id based on the subsplit bitset.
  size_t GetDAGNodeId(const Bitset &subsplit) const;
  // Gets the node id of the DAG root.
  size_t DAGRootNodeId() const;
  // Get the node based on the nodes id.
  SubsplitDAGNode *GetNode(const size_t node_id) const;
  // Return the node ids corresponding to the rootsplits.
  const SizeVector &GetRootsplitIds() const;
  // Get the GPCSP/edge index by its parent/child pair.
  size_t GetEdgeIdx(const Bitset &parent_subsplit, const Bitset &child_subsplit) const;
  size_t GetEdgeIdx(const size_t parent_id, const size_t child_id) const;
  // Get a sorted vector of all node subsplit bitset representation. Optionally only
  // graft nodes, or graft and host nodes.
  BitsetVector GetSortedVectorOfNodeBitsets(bool include_host = true) const;
  // Get a sorted vector of all edge PCSP bitset representation. Optionally only graft
  // edges, or graft and host edges.
  BitsetVector GetSortedVectorOfEdgeBitsets(bool include_host = true) const;

  // ** Counts

  // Total number of taxa in DAG (same as Subsplit length).
  size_t TaxonCount() const;
  // Total number of nodes in full proposed DAG.
  size_t NodeCount() const;
  // Total number of nodes in graft only.
  size_t GraftNodeCount() const;
  // Total number of nodes in host DAG only.
  size_t HostNodeCount() const;
  // Total number of edges in full proposed DAG.
  size_t EdgeCount() const;
  // Total number of edges in graft only.
  size_t GraftEdgeCount() const;
  // Total number of edges in host DAG only.
  size_t HostEdgeCount() const;
  // Check if a node is from the host, otherwise from the graft.
  bool IsNodeFromHost(size_t node_id) const;
  // Check if an edge is from the host, otherwise from the graft.
  bool IsEdgeFromHost(size_t edge_id) const;

  // ** Contains

  // Checks whether the node is in the host DAG or the graft.
  bool ContainsNode(const Bitset node_subsplit) const;
  bool ContainsNode(const size_t node_id) const;
  // Checks whether the node is in the graft only.
  bool ContainsGraftNode(const Bitset node_subsplit) const;
  bool ContainsGraftNode(const size_t node_id) const;
  // Checks whether the edge is in the host DAG or the graft.
  bool ContainsEdge(const size_t parent_id, const size_t child_id) const;
  // Checks whether the edge is in the graft only.
  bool ContainsGraftEdge(const size_t parent_id, const size_t child_id) const;

  // ** Validity Tests

  // Checks if host and graft are in a valid state.
  // Specifically, checks that:
  // - Host DAG is in valid state.
  // - Each node Id corresponds to its array index.
  // - Each node is valid. That is, either:
  //    - Node has zero parents and zero children.
  //    - Node has 1 or more parents, 1 or more right clade children, and 1 or more left
  //    clade children.
  // - All GraftDAG node ids are >= HostNodeCount.
  bool IsValid() const;
  // Checks if host and graft are in a valid state to add given node pair.
  // Specifically, check that:
  // - The nodes are adjacent.
  // - The nodes do not add/remove any taxa.
  // - The parent node has at least one parent.
  // - Including the child node, each clade of the parent node has at least one child.
  // - Each clade of the child node has at least 1 child.
  bool IsValidAddGraftNodePair(const Bitset parent_subsplit,
                               const Bitset child_subsplit) const;

 protected:
  // ** Modify DAG
  // These modifications do not ensure a valid, consistent state for DAG.

  // Add node to the graft.
  void CreateGraftNode(const Bitset &node_subsplit);
  // Add edge to the graft.
  void CreateGraftEdge(const size_t parent_id, const size_t child_id);
  // Add edge to the graft (Overload for when edge relation is known).
  void CreateGraftEdge(const size_t parent_id, const size_t child_id,
                       const bool is_left);
  // Connect main node to all adjacent nodes in vector.  `ignored_node_id_opt` provides
  // the option for one adjacent node that will not be connected.
  void ConnectNodeToAdjacentHostNodes(
      const size_t main_node_id, const SizeVector adjacent_node_ids,
      const bool is_main_node_parent, const bool is_left,
      std::optional<size_t> ignored_node_id_opt = std::nullopt);
  // Connect node in graft to all viable nodes in host.
  void ConnectNodeToHost(const size_t node_id, const Bitset node_subsplit,
                         std::optional<size_t> node_to_ignore = std::nullopt);

 protected:
  // DAG that the graft is proposed to be connected to.
  SubsplitDAG &host_dag_;
  // Nodes and edges in the graft.
  SubsplitDAGStorage graft_storage_;
  // Map of all DAG Nodes:
  //   - [ Node Subsplit (Bitset) => Node Id ]
  BitsetSizeMap subsplit_to_id_;
};
