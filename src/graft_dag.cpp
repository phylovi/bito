// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "graft_dag.hpp"

#include "gp_dag.hpp"
#include "subsplit_dag.hpp"

// ** Constructors

GraftDAG::GraftDAG(SubsplitDAG &dag) : host_dag_(dag) {}

// ** Comparators

int GraftDAG::Compare(const GraftDAG &other) const {
  return GraftDAG::Compare(*this, other);
}

int GraftDAG::Compare(const GraftDAG &lhs, const GraftDAG &rhs) {
  const SubsplitDAG &lhs_host = lhs.GetHostDAG();
  const SubsplitDAG &rhs_host = rhs.GetHostDAG();
  // (1) Compare host DAGs.
  int host_diff = SubsplitDAG::Compare(lhs_host, rhs_host);
  if (host_diff != 0) {
    return host_diff;
  }
  // (2) Compare graft nodes.
  const BitsetVector &lhs_nodes = lhs.GetSortedVectorOfNodeBitsets();
  const BitsetVector &rhs_nodes = rhs.GetSortedVectorOfNodeBitsets();
  if (lhs_nodes != rhs_nodes) {
    return (lhs_nodes < rhs_nodes) ? -1 : 1;
  }
  // (3) Compare graft edges.
  const BitsetVector &lhs_edges = lhs.GetSortedVectorOfEdgeBitsets();
  const BitsetVector &rhs_edges = rhs.GetSortedVectorOfEdgeBitsets();
  if (lhs_edges != rhs_edges) {
    return (lhs_edges < rhs_edges) ? -1 : 1;
  }
  return 0;
}

int GraftDAG::CompareToDAG(const SubsplitDAG &other) const {
  return GraftDAG::CompareToDAG(*this, other);
}

int GraftDAG::CompareToDAG(const GraftDAG &lhs, const SubsplitDAG &rhs) {
  // Compare taxon counts.
  const int taxon_diff = lhs.TaxonCount() - rhs.TaxonCount();
  if (taxon_diff != 0) {
    return taxon_diff;
  }
  // Compare nodes.
  const BitsetVector &lhs_nodes = lhs.GetSortedVectorOfNodeBitsets();
  const BitsetVector &rhs_nodes = rhs.GetSortedVectorOfNodeBitsets();
  if (lhs_nodes != rhs_nodes) {
    return (lhs_nodes < rhs_nodes) ? -1 : 1;
  }
  // Compare edges.
  const BitsetVector &lhs_edges = lhs.GetSortedVectorOfEdgeBitsets();
  const BitsetVector &rhs_edges = rhs.GetSortedVectorOfEdgeBitsets();
  if (lhs_edges != rhs_edges) {
    return (lhs_edges < rhs_edges) ? -1 : 1;
  }
  return 0;
}

// ** Modify GraftDAG

void GraftDAG::CreateGraftNode(const Bitset &node_subsplit) {
  // Do not add node to DAG if already exists in DAG.
  if (ContainsNode(node_subsplit)) {
    return;
  }
  // else, add node to graft node list.
  const size_t node_id = NodeCount();
  graft_storage_.AddVertex({node_id, node_subsplit}, HostNodeCount());
  // insert node into subsplit map.
  SafeInsert(subsplit_to_id_, node_subsplit, node_id);
}

void GraftDAG::CreateGraftEdge(const size_t parent_id, const size_t child_id) {
  // find relationship between nodes.
  const Bitset &parent_subsplit = GetDAGNode(parent_id).GetBitset();
  const Bitset &child_subsplit = GetDAGNode(child_id).GetBitset();
  bool is_left =
      !Bitset::SubsplitIsChildOfWhichParentClade(parent_subsplit, child_subsplit);
  CreateGraftEdge(parent_id, child_id, is_left);
}

void GraftDAG::CreateGraftEdge(const size_t parent_id, const size_t child_id,
                               const bool is_left) {
  // Add edge to edge list.
  Assert(ContainsNode(parent_id), "Node with the given parent_id does not exist.");
  Assert(ContainsNode(child_id), "Node with the given child_id does not exist.");

  auto parent_node = GetDAGNode(parent_id);
  auto child_node = GetDAGNode(child_id);

  // Add edge to grafted edges.
  if (!ContainsEdge(parent_id, child_id)) {
    const size_t &edge_id = EdgeCount();
    graft_storage_.AddLine(
        {edge_id, parent_id, child_id, is_left ? Clade::Left : Clade::Right},
        edge_id - graft_storage_.GetLines().size());
    // Add edge to individual SubsplitDAGNodes (only adds them to graft nodes).
    const bool parent_node_is_in_dag = !IsNodeFromHost(parent_id);
    const bool child_node_is_in_dag = !IsNodeFromHost(child_id);
    if (is_left) {
      if (!parent_node_is_in_dag) {
        parent_node.AddLeftLeafward(child_node.Id(), edge_id);
      }
      if (!child_node_is_in_dag) {
        child_node.AddLeftRootward(parent_node.Id(), edge_id);
      }
    } else {
      if (!parent_node_is_in_dag) {
        parent_node.AddRightLeafward(child_node.Id(), edge_id);
      }
      if (!child_node_is_in_dag) {
        child_node.AddRightRootward(parent_node.Id(), edge_id);
      }
    }
  }
}

void GraftDAG::AddGraftNodePair(const NNIOperation &nni) {
  return AddGraftNodePair(nni.parent_, nni.child_);
}

void GraftDAG::AddGraftNodePair(const Bitset &parent_subsplit,
                                const Bitset &child_subsplit) {
  Assert(IsValidAddGraftNodePair(parent_subsplit, child_subsplit),
         "AddGraftNodePair(): Is not valid node pair.");
  const bool parent_is_new =
      !(host_dag_.ContainsNode(parent_subsplit) || ContainsGraftNode(parent_subsplit));
  const bool child_is_new =
      !(host_dag_.ContainsNode(child_subsplit) || ContainsGraftNode(child_subsplit));
  // Assert that parent and child don't already exist in the DAG.
  if (!parent_is_new && !child_is_new) {
    return;
  }
  // Find relationship between nodes.
  const bool is_left = child_subsplit.SubsplitIsLeftChildOf(parent_subsplit);
  // Add nodes.
  if (parent_is_new) {
    CreateGraftNode(parent_subsplit);
  }
  const size_t parent_id = GetDAGNodeId(parent_subsplit);
  if (child_is_new) {
    CreateGraftNode(child_subsplit);
  }
  const size_t child_id = GetDAGNodeId(child_subsplit);
  // Find all adjacent host nodes and add edges to fully connect to new graft nodes.
  if (parent_is_new) {
    ConnectNodeToHost(parent_id, parent_subsplit, child_id);
  }
  if (child_is_new) {
    ConnectNodeToHost(child_id, child_subsplit, parent_id);
  }
  // Connect parent to child.
  CreateGraftEdge(parent_id, child_id, is_left);
}

void GraftDAG::RemoveAllGrafts() {
  graft_storage_ = {};
  subsplit_to_id_.clear();
}

void GraftDAG::ConnectNodeToHost(const size_t node_id, const Bitset node_subsplit,
                                 std::optional<size_t> node_to_ignore) {
  const auto [left_parents, right_parents] =
      host_dag_.BuildParentIdVectors(node_subsplit);
  ConnectNodeToAdjacentHostNodes(node_id, left_parents, false, true, node_to_ignore);
  ConnectNodeToAdjacentHostNodes(node_id, right_parents, false, false, node_to_ignore);
  const auto [left_children, right_children] =
      host_dag_.BuildChildIdVectors(node_subsplit);
  ConnectNodeToAdjacentHostNodes(node_id, left_children, true, true, node_to_ignore);
  ConnectNodeToAdjacentHostNodes(node_id, right_children, true, false, node_to_ignore);
}

void GraftDAG::ConnectNodeToAdjacentHostNodes(
    const size_t main_node_id, const SizeVector adjacent_node_ids,
    const bool is_main_node_parent, const bool is_left,
    std::optional<size_t> ignored_node_id_opt) {
  for (const auto adjacent_node_id : adjacent_node_ids) {
    if (ignored_node_id_opt && (adjacent_node_id == *ignored_node_id_opt)) {
      continue;
    }
    const size_t parent_id = (is_main_node_parent ? main_node_id : adjacent_node_id);
    const size_t child_id = (is_main_node_parent ? adjacent_node_id : main_node_id);
    CreateGraftEdge(parent_id, child_id, is_left);
  }
}

// ** Build Indexers/Vectors

std::pair<SizeVector, SizeVector> GraftDAG::BuildParentIdVectors(
    const Bitset &subsplit) const {
  SizeVector left_parents, right_parents;
  const auto &subsplit_map = host_dag_.GetSubsplitToIdMap();
  // Linear search for all parents.
  for (const auto &[potential_parent_subsplit, node_id] : subsplit_map) {
    if (subsplit.SubsplitIsLeftChildOf(potential_parent_subsplit)) {
      left_parents.push_back(node_id);
    } else if (subsplit.SubsplitIsRightChildOf(potential_parent_subsplit)) {
      right_parents.push_back(node_id);
    }
  }
  return {left_parents, right_parents};
}

std::pair<SizeVector, SizeVector> GraftDAG::BuildChildIdVectors(
    const Bitset &subsplit) const {
  SizeVector left_children, right_children;
  const auto &subsplit_map = host_dag_.GetSubsplitToIdMap();
  // Linear search for all parents.
  for (const auto &[potential_child_subsplit, node_id] : subsplit_map) {
    if (potential_child_subsplit.SubsplitIsLeftChildOf(subsplit)) {
      left_children.push_back(node_id);
    } else if (potential_child_subsplit.SubsplitIsRightChildOf(subsplit)) {
      right_children.push_back(node_id);
    }
  }
  return {left_children, right_children};
}

// ** Getters

const SubsplitDAG &GraftDAG::GetHostDAG() const { return host_dag_; }

SubsplitDAGNode GraftDAG::GetDAGNode(const size_t node_id) const {
  Assert(node_id < NodeCount(), "Node Id is out of valid range.");
  // Check if node is in the host DAG.
  if (node_id < HostNodeCount()) {
    return host_dag_.GetDAGNode(node_id);
  }
  // else, node is in the graft.
  return graft_storage_.GetVertices().at(node_id - host_dag_.NodeCount());
}

MutableSubsplitDAGNode GraftDAG::GetDAGNode(const size_t node_id) {
  Assert(node_id < NodeCount(), "Node Id is out of valid range.");
  // Check if node is in the host DAG.
  if (node_id < HostNodeCount()) {
    return host_dag_.GetDAGNode(node_id);
  }
  // else, node is in the graft.
  return graft_storage_.GetVertices().at(node_id - host_dag_.NodeCount());
}

size_t GraftDAG::GetDAGNodeId(const Bitset &node_subsplit) const {
  // Check if subsplit is in the host DAG.
  if (host_dag_.ContainsNode(node_subsplit)) {
    return host_dag_.GetDAGNodeId(node_subsplit);
  }
  // Check if subsplit is in the graft.
  Assert(subsplit_to_id_.find(node_subsplit) != subsplit_to_id_.end(),
         "GetDAGNodeId(): Subsplit does not correspond to a node in GraftDAG.");
  return subsplit_to_id_.at(node_subsplit);
}

size_t GraftDAG::GetDAGRootNodeId() const { return host_dag_.GetDAGRootNodeId(); }

size_t GraftDAG::GetEdgeIdx(const Bitset &parent_subsplit,
                            const Bitset &child_subsplit) const {
  const size_t &parent_id = GetDAGNodeId(parent_subsplit);
  const size_t &child_id = GetDAGNodeId(child_subsplit);
  size_t edge_idx = 0;
  if (host_dag_.ContainsEdge(parent_id, child_id)) {
    host_dag_.GetEdgeIdx(parent_id, child_id);
  } else {
    Assert(ContainsGraftEdge(parent_id, child_id),
           "GraftDAG::GetEdgeIdx(): Edge with given node pair does not exist.");
    edge_idx = graft_storage_.GetLine(parent_id, child_id).value().GetId();
  }
  return edge_idx;
}

size_t GraftDAG::GetEdgeIdx(const size_t parent_id, const size_t child_id) const {
  // Check for edge in graft.
  const auto edge = graft_storage_.GetLine(parent_id, child_id);
  if (edge.has_value()) {
    return edge.value().GetId();
  }
  // Check for edge in host.
  Assert(host_dag_.ContainsEdge(parent_id, child_id),
         "GetEdgeIdx: There is no edge corresponding to given parent/child pair in "
         "GraftDAG.");
  return host_dag_.GetEdgeIdx(parent_id, child_id);
}

BitsetVector GraftDAG::GetSortedVectorOfNodeBitsets(bool include_host) const {
  BitsetVector nodes;
  // Get graft node bitsets.
  for (const auto &subsplit_id_pair : subsplit_to_id_) {
    nodes.push_back(subsplit_id_pair.first);
  }
  // Get host node bitsets.
  if (include_host) {
    BitsetVector host_nodes = host_dag_.GetSortedVectorOfNodeBitsets();
    nodes.insert(nodes.end(), host_nodes.begin(), host_nodes.end());
  }
  std::sort(nodes.begin(), nodes.end());
  return nodes;
}

BitsetVector GraftDAG::GetSortedVectorOfEdgeBitsets(bool include_host) const {
  BitsetVector edges;
  // Get graft edge bitsets.
  for (const auto &[idpair_idx_pair, edge_idx] : graft_storage_.GetLines()) {
    std::ignore = edge_idx;
    const auto &[parent_id, child_id] = idpair_idx_pair;
    const auto &parent_bitset = GetDAGNode(parent_id).GetBitset();
    const auto &child_bitset = GetDAGNode(child_id).GetBitset();
    Bitset edge_bitset = Bitset::PCSP(parent_bitset, child_bitset);
    edges.push_back(edge_bitset);
  }
  // Get host edge bitsets.
  if (include_host) {
    BitsetVector host_edges = host_dag_.GetSortedVectorOfEdgeBitsets();
    edges.insert(edges.end(), host_edges.begin(), host_edges.end());
  }
  std::sort(edges.begin(), edges.end());

  return edges;
}

// ** Counts

size_t GraftDAG::TaxonCount() const { return host_dag_.TaxonCount(); }

size_t GraftDAG::NodeCount() const { return HostNodeCount() + GraftNodeCount(); }

size_t GraftDAG::GraftNodeCount() const { return graft_storage_.GetVertices().size(); }

size_t GraftDAG::HostNodeCount() const { return host_dag_.NodeCount(); }

size_t GraftDAG::EdgeCount() const { return HostEdgeCount() + GraftEdgeCount(); }

size_t GraftDAG::GraftEdgeCount() const { return graft_storage_.GetLines().size(); }

size_t GraftDAG::HostEdgeCount() const {
  return host_dag_.EdgeCountWithLeafSubsplits();
}

bool GraftDAG::IsNodeFromHost(size_t node_id) const {
  return (node_id < HostNodeCount());
}

bool GraftDAG::IsEdgeFromHost(size_t edge_id) const {
  return (edge_id < HostEdgeCount());
}

// ** Contains

bool GraftDAG::ContainsGraftNode(const Bitset node_subsplit) const {
  return subsplit_to_id_.find(node_subsplit) != subsplit_to_id_.end();
}

bool GraftDAG::ContainsGraftNode(const size_t node_id) const {
  return (node_id >= HostNodeCount()) && (node_id < NodeCount());
}

bool GraftDAG::ContainsNode(const Bitset node_subsplit) const {
  return host_dag_.ContainsNode(node_subsplit) || ContainsGraftNode(node_subsplit);
}

bool GraftDAG::ContainsNode(const size_t node_id) const {
  return host_dag_.ContainsNode(node_id) || ContainsGraftNode(node_id);
}

bool GraftDAG::ContainsGraftEdge(const size_t parent_id, const size_t child_id) const {
  return graft_storage_.GetLine(parent_id, child_id, HostNodeCount(), HostEdgeCount())
      .has_value();
}

bool GraftDAG::ContainsEdge(const size_t parent_id, const size_t child_id) const {
  return host_dag_.ContainsEdge(parent_id, child_id) ||
         ContainsGraftEdge(parent_id, child_id);
}

// ** Validation Tests

bool GraftDAG::IsValid() const {
  if (!host_dag_.IsValid()) {
    return false;
  }
  for (size_t i = 0; i < GraftNodeCount(); i++) {
    const size_t correct_id = HostNodeCount() + i;
    const auto node = graft_storage_.GetVertices()[i];
    if (node.Id() != correct_id) {
      return false;
    }
    if (!node.IsValid()) {
      return false;
    }
  }
  return true;
}

bool GraftDAG::IsValidAddGraftNodePair(const Bitset parent_subsplit,
                                       const Bitset child_subsplit) const {
  // Assert node pair does not already exist in host or graft.
  if (ContainsNode(parent_subsplit) && ContainsNode(child_subsplit)) {
    return false;
  }
  return host_dag_.IsValidAddNodePair(parent_subsplit, child_subsplit);
}
