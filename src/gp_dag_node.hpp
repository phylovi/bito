// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A class for node in a directed acyclic graph for likelihood EM training.

#ifndef SRC_GP_DAG_NODE_HPP_
#define SRC_GP_DAG_NODE_HPP_

#include <string>
#include <vector>

#include "bitset.hpp"

class GPDAGNode {
 public:
  GPDAGNode(size_t id, const Bitset &subsplit) : id_(id), subsplit_(subsplit) {}

  size_t Id() const { return id_; }
  const Bitset &GetBitset() const { return subsplit_; }
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
  const std::vector<size_t> &GetLeafwardRotated() { return leafward_rotated; }
  const std::vector<size_t> &GetLeafwardSorted() { return leafward_sorted; }
  const std::vector<size_t> &GetRootwardRotated() { return rootward_rotated; }
  const std::vector<size_t> &GetRootwardSorted() { return rootward_sorted; }

  std::string ToString();

 private:
  size_t id_;
  Bitset subsplit_;

  std::vector<size_t> leafward_rotated;
  std::vector<size_t> leafward_sorted;
  std::vector<size_t> rootward_rotated;
  std::vector<size_t> rootward_sorted;
};

#endif  // SRC_GP_DAG_NODE_HPP_
