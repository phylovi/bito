// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A node in a directed acyclic graph representing a collection of subsplits with their
// corresponding parent-child relationships.
//
// Each node represents a sorted osubsplit, which is stored as a bitset in `subsplit_`.
// The leafward edges are divided into two groups based on if they split apart the
// right clade (in which case they are called sorted children) or if they split apart
// the left clade (in which case they are called rotated children).
// Similarly, the rootward edges are divided into two groups based on if the child of
// the edge splits apart the right clade of the parent (in which case they are called
// sorted parents) or if they split apart the left clade of the parent (in which case
// they are called rotated parents).

#pragma once

#include "bitset.hpp"
#include "reindexer.hpp"
#include "sugar.hpp"

class SubsplitDAGNode {
 public:
  enum class Direction : bool {
    Rootward = 0,
    Leafward = 1,
  };
  enum class ParentClade : bool {
    Left = 0, 
    Right = 1, 
  };

  SubsplitDAGNode(size_t id, Bitset subsplit)
      : id_(id), subsplit_(std::move(subsplit)) {}

  size_t Id() const { return id_; }
  const Bitset &GetBitset() const { return subsplit_; }
  const Bitset GetBitset(bool rotated) const {
    return rotated ? subsplit_.SubsplitRotate() : subsplit_;
  }

  // ** Node Types

  // Is the node the DAG root (universal ancestor of the DAG)?
  bool IsDAGRootNode() const {
    return (right_rootward_.empty() && left_rootward_.empty());
  }
  // Is the node a rootsplit (direct descendent of root, dividing entire taxon set)?
  bool IsRootsplit() const { return subsplit_.SubsplitIsRootsplit(); }
  // Is the node a leaf (has no descendents)?
  bool IsLeaf() const {
    return left_leafward_.empty() && right_leafward_.empty();
  }

  // ** Edges 
  
  // Add edge from this node to adjacent_node.
  void AddEdge(size_t adjacent_node_id, bool is_leafward, bool is_rotated) {
    if (is_leafward) {
      is_rotated ? AddLeafwardLeftward(adjacent_node_id)
                 : AddLeafwardRightward(adjacent_node_id);
    } else {
      is_rotated ? AddRootwardLeftward(adjacent_node_id)
                 : AddRootwardRightward(adjacent_node_id);
    }
  }
  void AddEdge(size_t adjacent_node_id, Direction direction, ParentClade clade) {
    bool is_leafward = (direction == Direction::Leafward);
    bool is_rotated = (clade == ParentClade::Left);
    AddEdge(adjacent_node_id, is_leafward, is_rotated);
  }
  void AddLeafwardLeftward(size_t node_id) { left_leafward_.push_back(node_id); }
  void AddLeafwardRightward(size_t node_id) { right_leafward_.push_back(node_id); }
  void AddRootwardLeftward(size_t node_id) { left_rootward_.push_back(node_id); }
  void AddRootwardRightward(size_t node_id) { right_rootward_.push_back(node_id); }

  // Get vector of all adjacent node vectors along the specified direction.
  const SizeVector &GetEdge(bool is_leafward, bool is_rotated) {
    if (is_leafward) {
      auto &edges = is_rotated ? GetLeafwardLeftward() : GetLeafwardRightward();
      return edges;
    } else {
      auto &edges = is_rotated ? GetRootwardLeftward() : GetRootwardRightward();
      return edges;
    }
  }
  const SizeVector &GetEdge(Direction direction, ParentClade clade) {
    bool is_leafward = (direction == Direction::Leafward);
    bool is_rotated = (clade == ParentClade::Left);
    return GetEdge(is_leafward, is_rotated);
  }
  const SizeVector &GetLeafwardOrRootward(bool leafward, bool rotated) const {
    return leafward ? GetLeafward(rotated) : GetRootward(rotated);
  };
  const SizeVector &GetLeafwardLeftward() const { return left_leafward_; }
  const SizeVector &GetLeafwardRightward() const { return right_leafward_; }
  const SizeVector &GetLeafward(bool rotated) const {
    return rotated ? GetLeafwardLeftward() : GetLeafwardRightward();
  }
  const SizeVector &GetRootwardLeftward() const { return left_rootward_; }
  const SizeVector &GetRootwardRightward() const { return right_rootward_; }
  const SizeVector &GetRootward(bool rotated) const {
    return rotated ? GetRootwardLeftward() : GetRootwardRightward();
  }

  // After modifying parent DAG.  
  void RemapNodeIds(const SizeVector node_reindexer) {
    id_ = node_reindexer.at(id_);
    Reindexer::RemapIdVector(left_leafward_, node_reindexer);
    Reindexer::RemapIdVector(right_leafward_, node_reindexer);
    Reindexer::RemapIdVector(left_rootward_, node_reindexer);
    Reindexer::RemapIdVector(right_rootward_, node_reindexer);
  }

  bool IsValid() const;
  std::string ToString() const;

 private:
  // Node unique identifier.  Corresponds to the positional index in SubsplitDAG's
  // dag_nodes_.
  size_t id_;
  // Node bitset subsplit clades.
  const Bitset subsplit_;

  // List of adjacent nodes in all directions.
  SizeVector left_leafward_;
  SizeVector right_leafward_;
  SizeVector left_rootward_;
  SizeVector right_rootward_;
};
