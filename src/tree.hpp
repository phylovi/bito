// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_HPP_
#define SRC_TREE_HPP_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "node.hpp"
#include "typedefs.hpp"

class Tree {
 public:
  typedef std::shared_ptr<Tree> TreePtr;
  typedef std::vector<TreePtr> TreePtrVector;

  explicit Tree(Node::NodePtr root, TagDoubleMap branch_lengths)
      : root_(root), branch_lengths_(branch_lengths) {}

  const Node::NodePtr Root() const { return root_; }
  const TagDoubleMap BranchLengths() const { return branch_lengths_; }
  uint32_t LeafCount() const { return Root()->LeafCount(); }
  Node::NodePtrVec Children() const { return Root()->Children(); }

  std::string Newick(
      TagStringMapOption node_labels = std::experimental::nullopt) const {
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

  static TreePtr UnitBranchLengthTreeOf(Node::NodePtr topology) {
    TagDoubleMap branch_lengths;
    topology->PreOrder([&branch_lengths](const Node* node) {
      assert(branch_lengths.insert({node->Tag(), 1.}).second);
    });
    return std::make_shared<Tree>(topology, branch_lengths);
  }

  static TreePtrVector ExampleTrees() {
    TreePtrVector v;
    for (const auto& topology : Node::ExampleTopologies()) {
      v.push_back(UnitBranchLengthTreeOf(topology));
    }
    return v;
  }

 private:
  Node::NodePtr root_;
  TagDoubleMap branch_lengths_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_TREE_HPP_
