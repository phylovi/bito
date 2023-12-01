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
  std::map<NodeId, EdgeId> remapped;
  for (auto id = neighbors.begin(); id != neighbors.end(); ++id) {
    remapped[NodeId(node_reindexer[(*id).value_])] = EdgeId(id.GetEdgeId());
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

  void AddEdge(NodeId adjacent_node_id, EdgeId edge_id, Direction which_direction,
               SubsplitClade which_clade) {
    node_.AddNeighbor(which_direction, which_clade, adjacent_node_id, edge_id);
  }
  void AddLeftLeafward(NodeId node_id, EdgeId edge_id) {
    node_.AddNeighbor(Direction::Leafward, SubsplitClade::Left, node_id, edge_id);
  }
  void AddRightLeafward(NodeId node_id, EdgeId edge_id) {
    node_.AddNeighbor(Direction::Leafward, SubsplitClade::Right, node_id, edge_id);
  }
  void AddLeftRootward(NodeId node_id, EdgeId edge_id) {
    node_.AddNeighbor(Direction::Rootward, SubsplitClade::Left, node_id, edge_id);
  }
  void AddRightRootward(NodeId node_id, EdgeId edge_id) {
    node_.AddNeighbor(Direction::Rootward, SubsplitClade::Right, node_id, edge_id);
  }

  // ** Neighbors

  // This is a super-iterator that iterates over all ConstNeighborView iterators of all
  // given directions/clade combinations.
  class MultiConstNeighborsView {
   public:
    MultiConstNeighborsView(
        const SubsplitDAGNode& node,
        std::vector<std::pair<Direction, SubsplitClade>> dir_clade_pairs)
        : node_id_(node.Id()), dir_clade_pairs_(dir_clade_pairs) {
      for (const auto [dir, clade] : dir_clade_pairs_) {
        views_.push_back(node.GetNeighbors(dir, clade));
      }
    }

    class Iterator {
     public:
      Iterator(size_t view_idx, const NodeId node_id,
               const std::vector<std::pair<Direction, SubsplitClade>>& dir_clade_pairs,
               const std::vector<ConstNeighborsView>& views)
          : node_id_(node_id),
            dir_clade_pairs_(dir_clade_pairs),
            views_(views),
            view_idx_(view_idx) {
        AssignIterator();
        if (it_ != nullptr) {
          MaybeGetNextIterator();
        }
      }

      bool operator!=(const Iterator& other) { return view_idx_ != other.GetViewIdx(); }

      Iterator& operator++() {
        ++(*it_);
        MaybeGetNextIterator();
        return *this;
      }

      NodeId operator*() { return *(*it_); }

      NodeId GetNodeId() const { return (*it_).GetNodeId(); }
      EdgeId GetEdgeId() const { return (*it_).GetEdgeId(); }
      Direction GetDirection() const { return dir_clade_pairs_[view_idx_].first; }
      SubsplitClade GetSubsplitClade() const {
        return dir_clade_pairs_[view_idx_].second;
      }
      const ConstNeighborsView GetCurrentView() const { return views_[view_idx_]; }
      size_t GetViewIdx() const { return view_idx_; }

      NodeId GetParentId() const {
        return (GetDirection() == Direction::Leafward) ? node_id_ : GetNodeId();
      };
      NodeId GetChildId() const {
        return (GetDirection() == Direction::Leafward) ? GetNodeId() : node_id_;
      };

     private:
      ConstNeighborsView::Iterator GetIterator() { return *it_; }
      void MaybeGetNextIterator() {
        while ((*it_) == (*end_) and view_idx_ < views_.size()) {
          view_idx_++;
          AssignIterator();
        }
      }
      void AssignIterator() {
        if (view_idx_ >= views_.size()) return;
        it_ = std::make_unique<ConstNeighborsView::Iterator>(GetCurrentView().begin());
        end_ = std::make_unique<ConstNeighborsView::Iterator>(GetCurrentView().end());
      }

      const NodeId node_id_;
      const std::vector<std::pair<Direction, SubsplitClade>>& dir_clade_pairs_;
      const std::vector<ConstNeighborsView>& views_;
      std::unique_ptr<ConstNeighborsView::Iterator> it_ = nullptr;
      std::unique_ptr<ConstNeighborsView::Iterator> end_ = nullptr;
      size_t view_idx_;
    };

    auto begin() const { return Iterator(0, node_id_, dir_clade_pairs_, views_); }
    auto end() const {
      return Iterator(dir_clade_pairs_.size(), node_id_, dir_clade_pairs_, views_);
    }
    size_t size() const {
      size_t size = 0;
      std::for_each(views_.begin(), views_.end(),
                    [&size](const auto view) { size += view.size(); });
      return size;
    }
    bool empty() const { return size() == 0; }

   private:
    const NodeId node_id_;
    const std::vector<std::pair<Direction, SubsplitClade>> dir_clade_pairs_;
    std::vector<ConstNeighborsView> views_;
  };

  // Get all neighbors.
  MultiConstNeighborsView GetNeighbors() const {
    std::vector<std::pair<Direction, SubsplitClade>> dir_clade_pairs;
    for (auto dir : DirectionEnum::Iterator()) {
      for (auto clade : SubsplitCladeEnum::Iterator()) {
        dir_clade_pairs.push_back({dir, clade});
      }
    }
    return {node_, dir_clade_pairs};
  }
  // Get all leafward/rootward neighbors.
  MultiConstNeighborsView GetNeighbors(Direction dir) const {
    std::vector<std::pair<Direction, SubsplitClade>> dir_clade_pairs;
    for (auto clade : SubsplitCladeEnum::Iterator()) {
      dir_clade_pairs.push_back({dir, clade});
    }
  }

  // Get neighbors in specified direction.
  ConstNeighborsView GetNeighbors(Direction direction, SubsplitClade clade) const {
    return node_.GetNeighbors(direction, clade);
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
  ConstNeighborsView GetLeafward(SubsplitClade clade) const {
    return (clade == SubsplitClade::Left) ? GetLeftLeafward() : GetRightLeafward();
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
  ConstNeighborsView GetRootward(SubsplitClade clade) const {
    return (clade == SubsplitClade::Left) ? GetLeftRootward() : GetRightRootward();
  }

  // Remap node ids after modifying DAG.
  void RemapNodeIds(const Reindexer& node_reindexer) {
    Assert(node_reindexer.IsValid(),
           "Reindexer must be valid in GenericSubsplitDAGNode::RemapNodeIds.");
    node_.SetId(NodeId(node_reindexer.GetNewIndexByOldIndex(node_.GetId().value_)));
    node_.GetNeighbors(Direction::Leafward, SubsplitClade::Left)
        .RemapNodeIds(node_reindexer);
    node_.GetNeighbors(Direction::Leafward, SubsplitClade::Right)
        .RemapNodeIds(node_reindexer);
    node_.GetNeighbors(Direction::Rootward, SubsplitClade::Left)
        .RemapNodeIds(node_reindexer);
    node_.GetNeighbors(Direction::Rootward, SubsplitClade::Right)
        .RemapNodeIds(node_reindexer);
  }
  // Remap edge ids after modifying DAG.
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

  std::optional<std::tuple<EdgeId, Direction, SubsplitClade>> FindNeighbor(NodeId id) {
    return node_.FindNeighbor(id);
  }

  void SetEdgeId(NodeId neighbor, EdgeId line) { node_.SetEdgeId(neighbor, line); }

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
    str += std::to_string((*i).value_) + " ";
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
  std::string str =
      std::to_string(Id().value_) + ": " + GetBitset().SubsplitToString() + "\n";
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
  storage.AddLine({EdgeId(0), NodeId(0), NodeId(1), SubsplitClade::Left});
  storage.AddLine({EdgeId(1), NodeId(0), NodeId(2), SubsplitClade::Right});
  storage.AddLine({EdgeId(2), NodeId(1), NodeId(3), SubsplitClade::Left});
  storage.AddLine({EdgeId(3), NodeId(1), NodeId(4), SubsplitClade::Right});

  storage.AddVertex(DAGVertex{}.SetId(NodeId(0)))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Left, NodeId(1), EdgeId(0))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Right, NodeId(2), EdgeId(1));
  storage.AddVertex(DAGVertex{}.SetId(NodeId(1)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Left, NodeId(0), EdgeId(0))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Left, NodeId(3), EdgeId(2))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Right, NodeId(4), EdgeId(3));
  storage.AddVertex(DAGVertex{}.SetId(NodeId(2)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Right, NodeId(0), EdgeId(1));
  storage.AddVertex(DAGVertex{}.SetId(NodeId(3)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Left, NodeId(1), EdgeId(2));
  storage.AddVertex(DAGVertex{}.SetId(NodeId(4)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Right, NodeId(1), EdgeId(3));
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
               .GetEdgeId(),
           2);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
