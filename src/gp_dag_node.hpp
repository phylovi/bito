// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A class for node in a directed acyclic graph for generalized pruning.

#ifndef SRC_GP_DAG_NODE_HPP_
#define SRC_GP_DAG_NODE_HPP_

#include <string>
#include <vector>

#include "bitset.hpp"
#include "sugar.hpp"

class GPDAGNode {
 public:
  GPDAGNode(size_t id, const Bitset &subsplit) : id_(id), subsplit_(subsplit) {}

  size_t Id() const { return id_; }
  const Bitset &GetBitset() const { return subsplit_; }
  const Bitset GetBitset(bool rotated) const {
    return rotated ? subsplit_.RotateSubsplit() : subsplit_;
  }
  bool IsRoot() const {
    return (rootward_sorted.size() + rootward_rotated.size()) == 0;
  };
  bool IsLeaf() const {
    return (leafward_rotated.size() == 0) && (leafward_sorted.size() == 0);
  }

  void AddLeafwardRotated(size_t node_id) { leafward_rotated.push_back(node_id); }
  void AddLeafwardSorted(size_t node_id) { leafward_sorted.push_back(node_id); }
  void AddRootwardRotated(size_t node_id) { rootward_rotated.push_back(node_id); }
  void AddRootwardSorted(size_t node_id) { rootward_sorted.push_back(node_id); }
  const SizeVector &GetLeafwardRotated() const { return leafward_rotated; }
  const SizeVector &GetLeafwardSorted() const { return leafward_sorted; }
  const SizeVector &GetRootwardRotated() const { return rootward_rotated; }
  const SizeVector &GetRootwardSorted() const { return rootward_sorted; }

  std::string ToString();

 private:
  size_t id_;
  Bitset subsplit_;

  SizeVector leafward_rotated;
  SizeVector leafward_sorted;
  SizeVector rootward_rotated;
  SizeVector rootward_sorted;
};

#endif  // SRC_GP_DAG_NODE_HPP_
