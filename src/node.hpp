// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// The Node class is how we express tree topologies.
//
// Nodes are immutable after construction except for the index. The index is
// provided for applications where it is useful to have the edges numbered with
// a contiguous set of integers, where the leaf edges have the smallest such
// integers. Because this integer assignment cannot be known as we are building
// up the tree, we must make a second reindexing pass through the tree, which
// must mutate state. However, this reindexing pass is itself deterministic, so
// doing it a second time will always give the same result.
//
// In summary, call Reindex after building your tree if you need to use the
// index. Note that Tree construction calls Reindex, if you are manually
// manipulating the topology make you do manipulations with that in mind.

#ifndef SRC_NODE_HPP_
#define SRC_NODE_HPP_

#include <limits.h>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "intpack.hpp"
#include "typedefs.hpp"

class Node {
 public:
  typedef std::shared_ptr<Node> NodePtr;
  typedef std::vector<NodePtr> NodePtrVec;
  typedef std::shared_ptr<NodePtrVec> NodePtrVecPtr;
  typedef std::unordered_map<NodePtr, uint32_t> TopologyCounter;

  // This is the type of functions that are used in the PCSS recursion
  // functions.
  //
  // The signature is in 4 parts, each of which describes the
  // position in the tree and then the direction. The 4 parts are the sister
  // clade, the focal clade, child 0, and child 1. False means down the tree
  // structure and true means up.
  // See `doc/pcss.svg` for a diagram of the PCSS traversal. In that file,
  // the first tree shows the terminology, and the subsequent trees show the
  // calls to f_root and f_internal.
  typedef std::function<void(const Node*, bool, const Node*, bool, const Node*,
                             bool, const Node*, bool)>
      PCSSFun;

 public:
  explicit Node(uint32_t leaf_id);
  explicit Node(NodePtrVec children, size_t index);

  size_t Index() const { return index_; }
  uint64_t Tag() const { return tag_; }
  std::string TagString() const { return StringOfPackedInt(this->tag_); }
  uint32_t MaxLeafID() const { return UnpackFirstInt(tag_); }
  uint32_t LeafCount() const { return UnpackSecondInt(tag_); }
  size_t Hash() const { return hash_; }
  bool IsLeaf() const { return children_.empty(); }
  NodePtrVec Children() const { return children_; }

  bool operator==(const Node& other);

  void PreOrder(std::function<void(const Node*)> f) const;

  // Iterate f through (parent, sister, node) for internal nodes using a
  // preorder traversal.
  void TriplePreOrderInternal(
      std::function<void(const Node*, const Node*, const Node*)> f) const;

  // Traversal for rooted pairs in an unrooted subtree in its traditional rooted
  // representation.
  // We take in two functions, f_root, and f_internal, each of which take three
  // edges.
  // We assume that f_root is symmetric in its last two arguments so that
  // f_root's signature actually looks like f_root(node0, {node1, node2}).
  // We apply f_root to the descendant edges like so: 012, 120, and 201. Because
  // f_root is symmetric in the last two arguments, we are going to get all of
  // the distinct calls of f.
  // At the internal nodes we cycle through triples of (parent, sister, node)
  // for f_internal.
  void TriplePreOrder(
      std::function<void(const Node*, const Node*, const Node*)> f_root,
      std::function<void(const Node*, const Node*, const Node*)> f_internal)
      const;

  // See the typedef of PCSSFun to understand the argument type to this
  // function.
  void PCSSPreOrder(PCSSFun f) const;
  void PostOrder(std::function<void(const Node*)> f) const;
  void LevelOrder(std::function<void(const Node*)> f) const;

  // These two functions take functions accepting triples of (node_index,
  // child0_index, child1_index) and apply them according to various traversals.
  void BinaryIndexPreOrder(const std::function<void(int, int, int)> f) const;
  void BinaryIndexPostOrder(const std::function<void(int, int, int)> f) const;

  // Iterate f through (parent, sister, node) for bifurcating trees using a
  // preorder traversal.
  void TriplePreOrderBifurcating(
      std::function<void(const Node*, const Node*, const Node*)> f) const;
  void TripleIndexPreOrderBifurcating(
      std::function<void(int, int, int)> f) const;

  // This function assigns indices to the nodes of the topology: the leaves get
  // their indices (which are contiguously numbered from 0 through the leaf
  // count -1) and the rest get ordered according to a postorder traversal. Thus
  // the root always has index equal to the number of nodes in the tree.
  //
  // This function returns a map that maps the tags to their indices.
  TagSizeMap Reindex();

  std::string Newick(
      const DoubleVectorOption& branch_lengths = std::experimental::nullopt,
      const TagStringMapOption& node_labels = std::experimental::nullopt,
      bool show_tags = false) const;

  std::string NewickAux(const DoubleVectorOption& branch_lengths,
                        const TagStringMapOption& node_labels,
                        bool show_tags) const;

  // ** Class methods
  static NodePtr Leaf(uint32_t id);
  static NodePtr Join(NodePtrVec children, size_t index = SIZE_MAX);
  static NodePtr Join(NodePtr left, NodePtr right, size_t index = SIZE_MAX);

  static NodePtrVec ExampleTopologies();

  // A "cryptographic" hash function from Stack Overflow (the std::hash function
  // appears to leave uint32_ts as they are, which doesn't work for our
  // application).
  // https://stackoverflow.com/a/12996028/467327
  static uint32_t SOHash(uint32_t x);

  // Bit rotation from Stack Overflow.
  // c is the amount by which we rotate.
  // https://stackoverflow.com/a/776523/467327
  static size_t SORotate(size_t n, uint32_t c);

 private:
  NodePtrVec children_;
  // See beginning of file for notes about the index.
  size_t index_;
  // The tag_ is a pair of packed integers representing (1) the maximum leaf ID
  // of the leaves below this node, and (2) the number of leaves below the node.
  uint64_t tag_;
  size_t hash_;

  // Make copy constructors private to eliminate copying.
  Node(const Node&);
  Node& operator=(const Node&);

  // This is a private PostOrder that can change the Node.
  void MutablePostOrder(std::function<void(Node*)> f);
};

// Compare NodePtrs by their Nodes.
inline bool operator==(const Node::NodePtr& lhs, const Node::NodePtr& rhs) {
  return *lhs == *rhs;
}

inline bool operator!=(const Node::NodePtr& lhs, const Node::NodePtr& rhs) {
  return !(lhs == rhs);
}

namespace std {
template <>
struct hash<Node::NodePtr> {
  size_t operator()(const Node::NodePtr& n) const { return n->Hash(); }
};
template <>
struct equal_to<Node::NodePtr> {
  bool operator()(const Node::NodePtr& lhs, const Node::NodePtr& rhs) const {
    return lhs == rhs;
  }
};
}  // namespace std

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Node header") {
  Node::NodePtrVec examples = Node::ExampleTopologies();
  Node::NodePtr t1 = examples[0];
  Node::NodePtr t1_twin = examples[1];
  Node::NodePtr t2 = examples[2];
  Node::NodePtr t3 = examples[3];
  // TODO(ematsen) add real test for TriplePreorder
  std::cout << "TriplePreOrder" << std::endl;
  std::cout << t2->Newick() << std::endl;
  auto print_triple = [](const Node* parent, const Node* sister,
                         const Node* node) {
    std::cout << parent->TagString() << ", " << sister->TagString() << ", "
              << node->TagString() << " " << std::endl;
  };
  t2->TriplePreOrder(print_triple, print_triple);

  // This is actually a non-trivial test (see note in Node constructor above),
  // which shows why we need bit rotation.
  CHECK_NE(t1->Hash(), t2->Hash());

  CHECK_EQ(t1, t1_twin);
  CHECK_NE(t1, t2);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_NODE_HPP_
