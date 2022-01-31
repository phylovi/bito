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
#include "subsplit_dag_storage.hpp"

class SubsplitDAGNode {
 public:
  SubsplitDAGNode(const DAGTraits::Vertex& storage) : storage_(storage) {}

  size_t Id() const { return storage_.GetId(); }
  const Bitset &GetBitset() const { return storage_.GetSubsplit(); }
  const Bitset GetBitset(bool rotated) const {
    return rotated ? GetBitset().SubsplitRotate() : GetBitset();
  }
  bool IsDAGRootNode() const {
    return (GetRootwardSorted().empty() && GetRootwardRotated().empty());
  }
  bool IsRootsplit() const { return GetBitset().SubsplitIsRootsplit(); }
  bool IsLeaf() const { return GetLeafwardRotated().empty() && GetLeafwardSorted().empty(); }

  // void AddLeafwardRotated(size_t node_id) { leafward_rotated_.push_back(node_id); }
  // void AddLeafwardSorted(size_t node_id) { leafward_sorted_.push_back(node_id); }
  // void AddRootwardRotated(size_t node_id) { rootward_rotated_.push_back(node_id); }
  // void AddRootwardSorted(size_t node_id) { rootward_sorted_.push_back(node_id); }
  // #350 use enumerated types for rotated?
  SizeVector GetLeafwardOrRootward(bool leafward, bool rotated) const {
    return leafward ? GetLeafward(rotated) : GetRootward(rotated);
  };
  SizeVector GetLeafwardRotated() const { return GetIds(Direction::Leafward, Clade::Left); }
  SizeVector GetLeafwardSorted() const { return GetIds(Direction::Leafward, Clade::Right); }
  SizeVector GetLeafward(bool rotated) const {
    return rotated ? GetLeafwardRotated() : GetLeafwardSorted();
  }
  SizeVector GetRootwardRotated() const { return GetIds(Direction::Rootward, Clade::Left); }
  SizeVector GetRootwardSorted() const { return GetIds(Direction::Rootward, Clade::Right); }
  SizeVector GetRootward(bool rotated) const {
    return rotated ? GetRootwardRotated() : GetRootwardSorted();
  }
  void RemapNodeIds(const SizeVector node_reindexer) {
    throw std::runtime_error("Not implemented");
    // TODO
    // id_ = node_reindexer.at(id_);
    // Reindexer::RemapIdVector(leafward_rotated_, node_reindexer);
    // Reindexer::RemapIdVector(leafward_sorted_, node_reindexer);
    // Reindexer::RemapIdVector(rootward_rotated_, node_reindexer);
    // Reindexer::RemapIdVector(rootward_sorted_, node_reindexer);
  }

  std::string ToString() const;

 private:

  SizeVector GetIds(Direction direction, Clade clade) const {
    SizeVector ids;
    for (auto& i : storage_.GetNeighbors().Get(direction, clade)) ids.push_back(i.first);
    return ids;
  }

  // size_t id_;
  // const Bitset subsplit_;
  const DAGTraits::Vertex& storage_;

  // SizeVector leafward_rotated_;
  // SizeVector leafward_sorted_;
  // SizeVector rootward_rotated_;
  // SizeVector rootward_sorted_;
};
