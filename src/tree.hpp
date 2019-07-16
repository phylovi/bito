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

  explicit Tree(Node::NodePtr topology, TagDoubleMap branch_lengths)
      : topology_(topology), branch_lengths_(branch_lengths) {}

  const Node::NodePtr Topology() const { return topology_; }
  const TagDoubleMap BranchLengths() const { return branch_lengths_; }
  uint32_t LeafCount() const { return Topology()->LeafCount(); }
  Node::NodePtrVec Children() const { return Topology()->Children(); }

  bool operator==(const Tree& other) {
    return (this->Topology() == other.Topology()) &&
           (this->BranchLengths() == other.BranchLengths());
  }

  std::string Newick(
      TagStringMapOption node_labels = std::experimental::nullopt) const {
    return Topology()->Newick(branch_lengths_, node_labels);
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

  // Remove trifurcation at the root and make it a bifurcation.
  // Given (s0:b0, s1:b1, s2:b3), we get (s0:b0, (s1:b1, s2:b2):0):0.
  TreePtr Detrifurcate() {
    if (Children().size() != 3) {
      std::cerr << "Detrifurcate given a non-trifurcating tree.\n";
      abort();
    }
    auto branch_lengths = BranchLengths();
    auto root12 = Node::Join(Children()[1], Children()[2]);
    assert(branch_lengths.insert({root12->Tag(), 0.}).second);
    auto rerooted_topology = Node::Join(Children()[0], root12);
    if (branch_lengths.find(rerooted_topology->Tag()) == branch_lengths.end()) {
      assert(branch_lengths.insert({rerooted_topology->Tag(), 0.}).second);
    }
    return std::make_shared<Tree>(rerooted_topology, branch_lengths);
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
  Node::NodePtr topology_;
  TagDoubleMap branch_lengths_;
};

// Compare TreePtrs by their Trees.
inline bool operator==(const Tree::TreePtr& lhs, const Tree::TreePtr& rhs) {
  return *lhs == *rhs;
}

inline bool operator!=(const Tree::TreePtr& lhs, const Tree::TreePtr& rhs) {
  return !(lhs == rhs);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Tree") {
  auto trees = Tree::ExampleTrees();
  CHECK_EQ(trees[0]->Detrifurcate()->Topology(), trees[3]->Topology());
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_TREE_HPP_
