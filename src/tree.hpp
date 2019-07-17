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
  typedef std::vector<double> BranchLengthVector;

  // Note: any missing branch lengths are set to zero.
  explicit Tree(Node::NodePtr topology, TagDoubleMap branch_lengths)
      : topology_(topology) {
    auto tag_index_map = topology->Reindex();
    branch_lengths_.resize(topology->Index());
    for (const auto& iter : tag_index_map) {
      auto tag = iter.first;
      auto index = iter.second;
      auto search = branch_lengths.find(tag);
      if (search != branch_lengths.end()) {
        assert(index < branch_lengths_.size());
        branch_lengths_[index] = search->second;
      } else {
        branch_lengths_[index] = 0.;
      }
    }
  }
  explicit Tree(Node::NodePtr topology, BranchLengthVector branch_lengths)
      : topology_(topology), branch_lengths_(branch_lengths) {
    assert(topology->Index() == branch_lengths.size());
  }

  const Node::NodePtr Topology() const { return topology_; }
  const BranchLengthVector BranchLengths() const { return branch_lengths_; }
  uint32_t LeafCount() const { return Topology()->LeafCount(); }
  Node::NodePtrVec Children() const { return Topology()->Children(); }
  size_t Index() const { return Topology()->Index(); }

  bool operator==(const Tree& other) {
    return (this->Topology() == other.Topology()) &&
           (this->BranchLengths() == other.BranchLengths());
  }

  std::string Newick(
      TagStringMapOption node_labels = std::experimental::nullopt) const {
    return Topology()->Newick(branch_lengths_, node_labels);
  }

  double BranchLength(const Node* node) const {
    assert(node->Index() < branch_lengths_.size());
    return branch_lengths_[node->Index()];
  }

  // Remove trifurcation at the root and make it a bifurcation.
  // Given (s0:b0, s1:b1, s2:b3):b4, we get (s0:b0, (s1:b1, s2:b2):0):0.
  // Note that we zero out the root branch length.
  TreePtr Detrifurcate() {
    if (Children().size() != 3) {
      std::cerr << "Detrifurcate given a non-trifurcating tree.\n";
      abort();
    }
    auto branch_lengths = BranchLengths();
    auto our_index = Index();
    auto root12 = Node::Join(Children()[1], Children()[2], our_index);
    branch_lengths[our_index] = 0.;
    auto rerooted_topology = Node::Join(Children()[0], root12, our_index + 1);
    branch_lengths.push_back(0.);
    return std::make_shared<Tree>(rerooted_topology, branch_lengths);
  }

  static TreePtr UnitBranchLengthTreeOf(Node::NodePtr topology) {
    topology->Reindex();
    BranchLengthVector branch_lengths(topology->Index());
    topology->PreOrder([&branch_lengths](const Node* node) {
      branch_lengths[node->Index()] = 1.;
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
  BranchLengthVector branch_lengths_;
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
