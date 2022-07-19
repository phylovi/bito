// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A node in a directed acyclic graph representing a collection of subsplits with their
// corresponding parent-child relationships.
//
// Each node represents a sorted subsplit, which is stored as a bitset in DAGVertex.
// The leafward edges are divided into two groups based on if they split apart the
// right clade (in which case they are called sorted children) or if they split apart
// the left clade (in which case they are called rotated children).
// Similarly, the rootward edges are divided into two groups based on if the child of
// the edge splits apart the right clade of the parent (in which case they are called
// sorted parents) or if they split apart the left clade of the parent (in which case
// they are called rotated parents). The return type [Const]NeighborsView of querying
// edges is a view object, that can be used similarly to a vector reference: can be
// iterated, and checked for size() and empty().

#pragma once

#include "bitset.hpp"
#include "reindexer.hpp"
#include "sugar.hpp"
#include "subsplit_dag_storage.hpp"

namespace {
template <typename T>
auto& GetStorage(const T& node) {
  return node.node_;
}
static inline void RemapNeighbors(NeighborsView neighbors,
                                  const SizeVector& node_reindexer) {
  std::map<VertexId, LineId> remapped;
  for (auto i = neighbors.begin(); i != neighbors.end(); ++i) {
    remapped[VertexId(node_reindexer[*i])] = LineId(i.GetEdge());
  }
  neighbors.SetNeighbors(remapped);
}
}  // namespace

template <typename T>
class GenericSubsplitDAGNode {
 public:
  GenericSubsplitDAGNode(T& node) : node_{node} {}
  GenericSubsplitDAGNode(const GenericSubsplitDAGNode<std::remove_const_t<T>>& other)
      : node_{other.node_} {}

  // Compare SubsplitDAGNodes by their ids.
  static int Compare(const SubsplitDAGNode& node_a, const SubsplitDAGNode& node_b);
  // Compare SubsplitDAGNodes by their subsplit representations.
  static int CompareBySubsplit(const SubsplitDAGNode& node_a,
                               const SubsplitDAGNode& node_b);

  NodeId Id() const { return NodeId(node_.GetId()); }
  const Bitset& GetBitset() const { return node_.GetSubsplit(); }
  const Bitset GetBitset(bool rotated) const {
    return rotated ? node_.GetSubsplit().SubsplitRotate() : node_.GetSubsplit();
  }

  // ** Node Types

  // Is the node the DAG root (universal ancestor of the DAG)?
  bool IsDAGRootNode() const {
    return (GetRightRootward().empty() && GetLeftRootward().empty());
  }
  // Is the node a rootsplit (direct descendent of root, dividing entire taxon set)?
  bool IsRootsplit() const { return node_.GetSubsplit().SubsplitIsRootsplit(); }
  // Is the node a leaf (has no descendants)?
  bool IsLeaf() const {
    return GetLeftLeafward().empty() && GetRightLeafward().empty();
  }

  // ** Edges

  void AddEdge(NodeId adjacent_node_id, bool is_leafward, bool is_left,
               LineId edge_id) {
    if (is_leafward) {
      is_left ? AddLeftLeafward(adjacent_node_id, edge_id)
              : AddRightLeafward(adjacent_node_id, edge_id);
    } else {
      is_left ? AddLeftRootward(adjacent_node_id, edge_id)
              : AddRightRootward(adjacent_node_id, edge_id);
    }
  }

  void AddEdge(NodeId adjacent_node_id, Direction which_direction,
               SubsplitClade which_clade, LineId edge_id) {
    bool is_leafward = (which_direction == Direction::Leafward);
    bool is_left = (which_clade == SubsplitClade::Left);
    AddEdge(adjacent_node_id, is_leafward, is_left, edge_id);
  }

  void AddLeftLeafward(NodeId node_id, LineId edge_id) {
    node_.AddNeighbor(Direction::Leafward, SubsplitClade::Left, node_id, edge_id);
  }
  void AddRightLeafward(NodeId node_id, LineId edge_id) {
    node_.AddNeighbor(Direction::Leafward, SubsplitClade::Right, node_id, edge_id);
  }
  void AddLeftRootward(NodeId node_id, LineId edge_id) {
    node_.AddNeighbor(Direction::Rootward, SubsplitClade::Left, node_id, edge_id);
  }
  void AddRightRootward(NodeId node_id, LineId edge_id) {
    node_.AddNeighbor(Direction::Rootward, SubsplitClade::Right, node_id, edge_id);
  }

  // Get a view into all adjacent node vectors along the specified direction.
  ConstNeighborsView GetEdge(bool is_leafward, bool is_left) const {
    if (is_leafward) {
      return is_left ? GetLeftLeafward() : GetRightLeafward();
    } else {
      return is_left ? GetLeftRootward() : GetRightRootward();
    }
  }

  ConstNeighborsView GetNeighbors(Direction direction, SubsplitClade clade) const {
    return node_.GetNeighbors(direction, clade);
  }

  void AddEdge(NodeId adjacent_node_id, LineId edge_id, Direction which_direction,
               SubsplitClade which_clade) {
    node_.AddNeighbor(which_direction, which_clade, adjacent_node_id, edge_id);
  }

  ConstNeighborsView GetLeafwardOrRootward(bool leafward, bool rotated) const {
    return leafward ? GetLeafward(rotated) : GetRootward(rotated);
  };
  ConstNeighborsView GetLeftLeafward() const {
    return node_.GetNeighbors(Direction::Leafward, SubsplitClade::Left);
  }
  ConstNeighborsView GetRightLeafward() const {
    return node_.GetNeighbors(Direction::Leafward, SubsplitClade::Right);
  }
  ConstNeighborsView GetLeafward(bool rotated) const {
    return rotated ? GetLeftLeafward() : GetRightLeafward();
  }
  ConstNeighborsView GetLeftRootward() const {
    return node_.GetNeighbors(Direction::Rootward, SubsplitClade::Left);
  }
  ConstNeighborsView GetRightRootward() const {
    return node_.GetNeighbors(Direction::Rootward, SubsplitClade::Right);
  }
  ConstNeighborsView GetRootward(bool rotated) const {
    return rotated ? GetLeftRootward() : GetRightRootward();
  }

  // Remap node ids modifying parent DAG.
  void RemapNodeIds(const Reindexer& node_reindexer) {
    Assert(node_reindexer.IsValid(),
           "Reindexer must be valid in GenericSubsplitDAGNode::RemapNodeIds.");
    node_.SetId(NodeId(node_reindexer.GetNewIndexByOldIndex(size_t(node_.GetId()))));
    node_.GetNeighbors(Direction::Leafward, SubsplitClade::Left)
        .RemapNodeIds(node_reindexer);
    node_.GetNeighbors(Direction::Leafward, SubsplitClade::Right)
        .RemapNodeIds(node_reindexer);
    node_.GetNeighbors(Direction::Rootward, SubsplitClade::Left)
        .RemapNodeIds(node_reindexer);
    node_.GetNeighbors(Direction::Rootward, SubsplitClade::Right)
        .RemapNodeIds(node_reindexer);
  }
  // Remap edge ids after modifying parent DAG.
  void RemapEdgeIdxs(const Reindexer& edge_reindexer) {
    Assert(edge_reindexer.IsValid(),
           "Reindexer must be valid in GenericSubsplitDAGNode::RemapNodeIds.");
    node_.GetNeighbors(Direction::Leafward, SubsplitClade::Left)
        .RemapEdgeIdxs(edge_reindexer);
    node_.GetNeighbors(Direction::Leafward, SubsplitClade::Right)
        .RemapEdgeIdxs(edge_reindexer);
    node_.GetNeighbors(Direction::Rootward, SubsplitClade::Left)
        .RemapEdgeIdxs(edge_reindexer);
    node_.GetNeighbors(Direction::Rootward, SubsplitClade::Right)
        .RemapEdgeIdxs(edge_reindexer);
  }

  std::optional<std::tuple<LineId, Direction, SubsplitClade>> FindNeighbor(
      VertexId id) {
    return node_.FindNeighbor(id);
  }

  void SetLineId(VertexId neighbor, LineId line) { node_.SetLineId(neighbor, line); }

  bool IsValid() const;
  std::string ToString() const;

 private:
  friend class DAGVertex;
  template <typename>
  friend class GenericSubsplitDAGNode;
  friend DAGVertex& GetStorage(const GenericSubsplitDAGNode<DAGVertex>&);
  friend const DAGVertex& GetStorage(const GenericSubsplitDAGNode<const DAGVertex>&);
  T& node_;
};

namespace {

static inline std::string GetNeighborString(ConstNeighborsView neighbors) {
  std::string str;
  for (auto i = neighbors.begin(); i != neighbors.end(); ++i) {
    str += std::to_string(*i) + " ";
  }
  return str;
}
}  // namespace

template <typename T>
bool GenericSubsplitDAGNode<T>::IsValid() const {
  // If node is a leaf, then a valid node should have no parents.
  if (IsLeaf()) {
    return (GetRightLeafward().size() + GetLeftLeafward().size() == 0);
  }
  // If node is a root, then a valid node should have no children.
  else if (IsDAGRootNode()) {
    return (GetRightRootward().size() + GetLeftRootward().size() == 0);
  }
  // If neither, then node should either have:
  // (1) Zero parents and zero children.
  // (2) 1+ parents, 1+ sorted children, and 1+ rotated children.
  size_t parent_node_count = GetRightRootward().size() + GetLeftRootward().size();
  if (parent_node_count > 0) {
    if (GetRightLeafward().size() == 0 || GetLeftLeafward().size() == 0) {
      return false;
    }
  } else {
    if (GetRightLeafward().size() > 0 || GetLeftLeafward().size() > 0) {
      return false;
    }
  }
  return true;
}

template <typename T>
std::string GenericSubsplitDAGNode<T>::ToString() const {
  std::string str = std::to_string(Id()) + ": " + GetBitset().SubsplitToString() + "\n";
  str += "Right Rootward: " + GetNeighborString(GetRightRootward()) + "\n";
  str += "Left Rootward: " + GetNeighborString(GetLeftRootward()) + "\n";
  str += "Right Leafward: " + GetNeighborString(GetRightLeafward()) + "\n";
  str += "Left Leafward: " + GetNeighborString(GetLeftLeafward()) + "\n";
  return str;
}

DAGVertex::DAGVertex(SubsplitDAGNode node) : DAGVertex(node.node_) {}
DAGVertex::DAGVertex(MutableSubsplitDAGNode node) : DAGVertex(node.node_) {}

template <typename T>
typename GenericVerticesView<T>::view_type GenericVerticesView<T>::Iterator::operator*()
    const {
  return storage_.vertices_[index_];
}
template <typename T>
typename GenericVerticesView<T>::view_type GenericVerticesView<T>::operator[](
    size_t i) const {
  return storage_.vertices_[i];
}
template <typename T>
typename GenericVerticesView<T>::view_type GenericVerticesView<T>::at(size_t i) const {
  return storage_.vertices_.at(i);
}

#ifdef DOCTEST_LIBRARY_INCLUDED

inline DAGVertex& GetStorage(const GenericSubsplitDAGNode<DAGVertex>& node) {
  return node.node_;
}
inline const DAGVertex& GetStorage(
    const GenericSubsplitDAGNode<const DAGVertex>& node) {
  return node.node_;
}

/* Create the following topology:
            [0]
            / \
           0   1
          /     \
        [1]     [2]
        / \
       2   3
      /     \
    [3]     [4]
 */

static SubsplitDAGStorage MakeStorage() {
  SubsplitDAGStorage storage;
  storage.AddLine({LineId(0), VertexId(0), VertexId(1), SubsplitClade::Left});
  storage.AddLine({LineId(1), VertexId(0), VertexId(2), SubsplitClade::Right});
  storage.AddLine({LineId(2), VertexId(1), VertexId(3), SubsplitClade::Left});
  storage.AddLine({LineId(3), VertexId(1), VertexId(4), SubsplitClade::Right});

  storage.AddVertex(DAGVertex{}.SetId(VertexId(0)))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Left, VertexId(1), LineId(0))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Right, VertexId(2), LineId(1));
  storage.AddVertex(DAGVertex{}.SetId(VertexId(1)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Left, VertexId(0), LineId(0))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Left, VertexId(3), LineId(2))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Right, VertexId(4), LineId(3));
  storage.AddVertex(DAGVertex{}.SetId(VertexId(2)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Right, VertexId(0), LineId(1));
  storage.AddVertex(DAGVertex{}.SetId(VertexId(3)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Left, VertexId(1), LineId(2));
  storage.AddVertex(DAGVertex{}.SetId(VertexId(4)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Right, VertexId(1), LineId(3));
  return storage;
}

TEST_CASE("SubsplitDAGStorage: LinesView structured binding") {
  auto storage = MakeStorage();

  size_t i = 0;
  for (auto [node_ids, line_id] : storage.GetLines()) {
    std::ignore = line_id;
    auto [parent_id, child_id] = node_ids;
    switch (i++) {
      case 0:
        CHECK_EQ(parent_id, 0);
        CHECK_EQ(child_id, 1);
        break;
      case 1:
        CHECK_EQ(parent_id, 0);
        CHECK_EQ(child_id, 2);
        break;
      case 2:
        CHECK_EQ(parent_id, 1);
        CHECK_EQ(child_id, 3);
        break;
      case 3:
        CHECK_EQ(parent_id, 1);
        CHECK_EQ(child_id, 4);
        break;
      default:
        Failwith("More lines than expected");
    }
  }
}

TEST_CASE("SubsplitDAGStorage: Neighbors iterator") {
  auto storage = MakeStorage();

  CHECK_EQ(*GetStorage(storage.GetVertices()[1])
                .GetNeighbors(Direction::Leafward, SubsplitClade::Left)
                .begin(),
           3);
  CHECK_EQ(GetStorage(storage.GetVertices()[1])
               .GetNeighbors(Direction::Leafward, SubsplitClade::Left)
               .begin()
               .GetEdge(),
           2);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
