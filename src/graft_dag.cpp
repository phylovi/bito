// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "graft_dag.hpp"

// ** Constructors

GraftDAG::GraftDAG(SubsplitDAG &dag)
    : SubsplitDAG{dag, HostDispatchTag{}}, host_dag_{dag} {}

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
  auto lhs_nodes = lhs.BuildSetOfNodeBitsets();
  auto rhs_nodes = rhs.BuildSetOfNodeBitsets();
  if (lhs_nodes != rhs_nodes) {
    return (lhs_nodes < rhs_nodes) ? -1 : 1;
  }
  // (3) Compare graft edges.
  auto lhs_edges = lhs.BuildSetOfEdgeBitsets();
  auto rhs_edges = rhs.BuildSetOfEdgeBitsets();
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
  auto lhs_nodes = lhs.BuildSetOfNodeBitsets();
  auto rhs_nodes = rhs.BuildSetOfNodeBitsets();
  if (lhs_nodes != rhs_nodes) {
    return (lhs_nodes < rhs_nodes) ? -1 : 1;
  }
  // Compare edges.
  auto lhs_edges = lhs.BuildSetOfEdgeBitsets();
  auto rhs_edges = rhs.BuildSetOfEdgeBitsets();
  if (lhs_edges != rhs_edges) {
    return (lhs_edges < rhs_edges) ? -1 : 1;
  }
  return 0;
}

// ** Modify GraftDAG

SubsplitDAG::ModificationResult GraftDAG::AddNodePair(const NNIOperation &nni) {
  return AddNodePair(nni.parent_, nni.child_);
}

SubsplitDAG::ModificationResult GraftDAG::AddNodePair(const Bitset &parent_subsplit,
                                                      const Bitset &child_subsplit) {
  // GetHostDAG().IsValidAddNodePair(parent_subsplit, child_subsplit);
  return AddNodes({{parent_subsplit, child_subsplit}});
}

SubsplitDAG::ModificationResult GraftDAG::AddNodes(
    const BitsetPairVector &node_subsplit_pairs) {
  auto mods = SubsplitDAG::AddNodePairInternals(node_subsplit_pairs);
  graft_node_count_ += mods.added_node_ids.size();
  graft_edge_count_ += mods.added_edge_idxs.size();
  return mods;
}

void GraftDAG::RemoveAllGrafts() {
  ResetHostDAG(host_dag_);
  graft_node_count_ = 0;
  graft_edge_count_ = 0;
}

// ** Getters

const SubsplitDAG &GraftDAG::GetHostDAG() const { return host_dag_; }

// ** Counts

size_t GraftDAG::GraftNodeCount() const {
  return storage_.GetVertices().size() - storage_.HostVerticesCount();
}

size_t GraftDAG::HostNodeCount() const { return storage_.HostVerticesCount(); }

size_t GraftDAG::GraftEdgeCount() const {
  return storage_.GetLines().size() - storage_.HostLinesCount();
}

size_t GraftDAG::HostEdgeCount() const {
  return GetHostDAG().EdgeCountWithLeafSubsplits();
}

bool GraftDAG::IsNodeFromHost(NodeId node_id) const {
  return node_id < HostNodeCount();
}

bool GraftDAG::IsEdgeFromHost(EdgeId edge_id) const {
  return edge_id < HostEdgeCount();
}

// ** Contains

bool GraftDAG::ContainsGraftNode(const Bitset node_subsplit) const {
  if (!ContainsNode(node_subsplit)) return false;
  return !IsNodeFromHost(GetDAGNodeId(node_subsplit));
}

bool GraftDAG::ContainsGraftNode(const NodeId node_id) const {
  if (IsNodeFromHost(node_id)) return false;
  return ContainsNode(node_id);
}

bool GraftDAG::ContainsGraftEdge(const NodeId parent_id, const NodeId child_id) const {
  auto edge = storage_.GetLine(parent_id, child_id);
  if (!edge.has_value()) return false;
  return !IsEdgeFromHost(edge.value().GetId());
}

bool GraftDAG::ContainsGraftEdge(const EdgeId edge_idx) const {
  auto edge = storage_.GetLine(edge_idx);
  if (!edge.has_value()) return false;
  return !IsEdgeFromHost(edge.value().GetId());
}

// ** Miscellaneous

size_t GraftDAG::GetPLVIndex(PLVNodeHandler::PLVType plv_type, NodeId node_id) const {
  return PLVNodeHandler::GetPVIndex(plv_type, node_id, NodeCountWithoutDAGRoot())
      .value_;
}
