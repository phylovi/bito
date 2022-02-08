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
    ROOTWARD = 0,
    LEAFWARD = 1,
  };
  enum class ParentClade : bool {
    LEFTSIDE = 0, /* ROTATED */
    RIGHTSIDE = 1 /* SORTED */
  };

  SubsplitDAGNode(size_t id, Bitset subsplit)
      : id_(id), subsplit_(std::move(subsplit)) {}

  size_t Id() const { return id_; }
  const Bitset &GetBitset() const { return subsplit_; }
  const Bitset GetBitset(bool rotated) const {
    return rotated ? subsplit_.SubsplitRotate() : subsplit_;
  }
  bool IsDAGRootNode() const {
    return (rootward_rightward_.empty() && rootward_leftward_.empty());
  }
  //
  bool IsRootsplit() const { return subsplit_.SubsplitIsRootsplit(); }
  bool IsLeaf() const {
    return leafward_leftward_.empty() && leafward_rightward_.empty();
  }

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
    bool is_leafward = (direction == Direction::LEAFWARD);
    bool is_rotated = (clade == ParentClade::LEFTSIDE);
    AddEdge(adjacent_node_id, is_leafward, is_rotated);
  }
  void AddLeafwardLeftward(size_t node_id) { leafward_leftward_.push_back(node_id); }
  void AddLeafwardRightward(size_t node_id) { leafward_rightward_.push_back(node_id); }
  void AddRootwardLeftward(size_t node_id) { rootward_leftward_.push_back(node_id); }
  void AddRootwardRightward(size_t node_id) { rootward_rightward_.push_back(node_id); }

  // #350 use enumerated types for rotated?
  // Get vector of all adjacent node vectors along the specified direction.
  void GetEdge(bool is_leafward, bool is_rotated) {
    if (is_leafward) {
      is_rotated ? GetLeafwardLeftward() : GetLeafwardRightward();
    } else {
      is_rotated ? GetRootwardLeftward() : GetRootwardRightward();
    }
  }
  void GetEdge(Direction direction, ParentClade clade) {
    bool is_leafward = (direction == Direction::LEAFWARD);
    bool is_rotated = (clade == ParentClade::LEFTSIDE);
    GetEdge(is_leafward, is_rotated);
  }
  const SizeVector &GetLeafwardOrRootward(bool leafward, bool rotated) const {
    return leafward ? GetLeafward(rotated) : GetRootward(rotated);
  };
  const SizeVector &GetLeafwardLeftward() const { return leafward_leftward_; }
  const SizeVector &GetLeafwardRightward() const { return leafward_rightward_; }
  const SizeVector &GetLeafward(bool rotated) const {
    return rotated ? GetLeafwardLeftward() : GetLeafwardRightward();
  }
  const SizeVector &GetRootwardLeftward() const { return rootward_leftward_; }
  const SizeVector &GetRootwardRightward() const { return rootward_rightward_; }
  const SizeVector &GetRootward(bool rotated) const {
    return rotated ? GetRootwardLeftward() : GetRootwardRightward();
  }
  void RemapNodeIds(const SizeVector node_reindexer) {
    id_ = node_reindexer.at(id_);
    Reindexer::RemapIdVector(leafward_leftward_, node_reindexer);
    Reindexer::RemapIdVector(leafward_rightward_, node_reindexer);
    Reindexer::RemapIdVector(rootward_leftward_, node_reindexer);
    Reindexer::RemapIdVector(rootward_rightward_, node_reindexer);
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
  SizeVector leafward_leftward_;
  SizeVector leafward_rightward_;
  SizeVector rootward_leftward_;
  SizeVector rootward_rightward_;
};
