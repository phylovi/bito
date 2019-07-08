// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_HPP_
#define SRC_TREE_HPP_

#include <memory>
#include <string>
#include <unordered_map>
#include "node.hpp"
#include "typedefs.hpp"

class Tree {
 public:
  typedef std::shared_ptr<Tree> TreePtr;

  explicit Tree(Node::NodePtr root, TagDoubleMap branch_lengths)
      : root_(root), branch_lengths_(branch_lengths) {}

  const Node::NodePtr Root() const { return root_; }
  uint32_t LeafCount() const { return Root()->LeafCount(); }
  std::string Newick(TagStringMapOption node_labels = std::nullopt) const {
    return Root()->Newick(branch_lengths_, node_labels);
  }
  double BranchLength(const Node* node) const {
    auto search = branch_lengths_.find(node->Tag());
    if (search != branch_lengths_.end()) {
      return search->second;
    } else {
      std::cerr << "Branch length not found for node tagged '"
                << node->TagString() << "'.\n";
      abort();
    }
  }

 private:
  Node::NodePtr root_;
  TagDoubleMap branch_lengths_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_HPP_
