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
}  // namespace

template <typename T>
class GenericSubsplitDAGNode {
 public:
  GenericSubsplitDAGNode(T& node) : node_{node} {}
  GenericSubsplitDAGNode(const GenericSubsplitDAGNode<std::remove_const_t<T>>& other)
      : node_{other.node_} {}

  size_t Id() const { return node_.GetId(); }
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

  // Add edge from this node to adjacent_node.
  void AddEdge(size_t adjacent_node_id, LineId edge_id, Direction which_direction,
               Clade which_clade) {
    node_.AddNeighbor(which_direction, which_clade, adjacent_node_id, edge_id);
  }

  ConstNeighborsView GetEdge(bool is_leafward, bool is_rotated) const {
    if (is_leafward) {
      return is_rotated ? GetLeftLeafward() : GetRightLeafward();
    } else {
      return is_rotated ? GetLeftRootward() : GetRightRootward();
    }
  }
  ConstNeighborsView GetEdge(Direction direction, Clade clade) const {
    bool is_leafward = (direction == Direction::Leafward);
    bool is_rotated = (clade == Clade::Left);
    return GetEdge(is_leafward, is_rotated);
  }
  ConstNeighborsView GetLeafwardOrRootward(bool leafward, bool rotated) const {
    return leafward ? GetLeafward(rotated) : GetRootward(rotated);
  };
  ConstNeighborsView GetLeftLeafward() const {
    return node_.GetNeighbors(Direction::Leafward, Clade::Left);
  }
  ConstNeighborsView GetRightLeafward() const {
    return node_.GetNeighbors(Direction::Leafward, Clade::Right);
  }
  ConstNeighborsView GetLeafward(bool rotated) const {
    return rotated ? GetLeftLeafward() : GetRightLeafward();
  }
  ConstNeighborsView GetLeftRootward() const {
    return node_.GetNeighbors(Direction::Rootward, Clade::Left);
  }
  ConstNeighborsView GetRightRootward() const {
    return node_.GetNeighbors(Direction::Rootward, Clade::Right);
  }
  ConstNeighborsView GetRootward(bool rotated) const {
    return rotated ? GetLeftRootward() : GetRightRootward();
  }

  // After modifying parent DAG.
  void RemapNodeIds(const SizeVector node_reindexer) {
    node_.SetId(node_reindexer.at(node_.GetId()));
    Reindexer::RemapIdVector(node_.GetNeighbors(Direction::Leafward, Clade::Left),
                             node_reindexer);
    Reindexer::RemapIdVector(node_.GetNeighbors(Direction::Leafward, Clade::Right),
                             node_reindexer);
    Reindexer::RemapIdVector(node_.GetNeighbors(Direction::Rootward, Clade::Left),
                             node_reindexer);
    Reindexer::RemapIdVector(node_.GetNeighbors(Direction::Rootward, Clade::Right),
                             node_reindexer);
  }

  std::optional<std::tuple<LineId, Direction, Clade>> FindNeighbor(VertexId id) {
    return node_.FindNeighbor(id);
  }

  void SetLineId(VertexId neighbor, LineId line) { node_.SetLineId(neighbor, line); }

  bool IsValid() const;
  std::string ToString() const;

 private:
  friend class DAGVertex;
  template <typename>
  friend class GenericSubsplitDAGNode;
  friend auto& GetStorage<>(const GenericSubsplitDAGNode&);
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
  storage.AddLine({0, 0, 1, Clade::Left});
  storage.AddLine({1, 0, 2, Clade::Right});
  storage.AddLine({2, 1, 3, Clade::Left});
  storage.AddLine({3, 1, 4, Clade::Right});

  storage.AddVertex(DAGVertex{}.SetId(0))
      .AddNeighbor(Direction::Leafward, Clade::Left, 1, 0)
      .AddNeighbor(Direction::Leafward, Clade::Right, 2, 1);
  storage.AddVertex(DAGVertex{}.SetId(1))
      .AddNeighbor(Direction::Rootward, Clade::Left, 0, 0)
      .AddNeighbor(Direction::Leafward, Clade::Left, 3, 2)
      .AddNeighbor(Direction::Leafward, Clade::Right, 4, 3);
  storage.AddVertex(DAGVertex{}.SetId(2))
      .AddNeighbor(Direction::Rootward, Clade::Right, 0, 1);
  storage.AddVertex(DAGVertex{}.SetId(3))
      .AddNeighbor(Direction::Rootward, Clade::Left, 1, 2);
  storage.AddVertex(DAGVertex{}.SetId(4))
      .AddNeighbor(Direction::Rootward, Clade::Right, 1, 3);
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
                .GetNeighbors(Direction::Leafward, Clade::Left)
                .begin(),
           3);
  CHECK_EQ(GetStorage(storage.GetVertices()[1])
               .GetNeighbors(Direction::Leafward, Clade::Left)
               .begin()
               .GetEdge(),
           2);

  *GetStorage(storage.GetVertices()[1])
       .GetNeighbors(Direction::Leafward, Clade::Left)
       .begin() = 42;

  CHECK_EQ(*GetStorage(storage.GetVertices()[1])
                .GetNeighbors(Direction::Leafward, Clade::Left)
                .begin(),
           42);
  CHECK_EQ(GetStorage(storage.GetVertices()[1])
               .GetNeighbors(Direction::Leafward, Clade::Left)
               .begin()
               .GetEdge(),
           2);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
