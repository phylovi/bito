// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// The Node class is how we express tree topologies.
//
// Nodes are immutable after construction except for the id_ and the leaves_.
// The id_ is provided for applications where it is useful to have the edges
// numbered with a contiguous set of integers. The leaves get
// their indices (which are contiguously numbered from 0 through the leaf
// count minus 1) and the rest get ordered according to a postorder traversal.
// Thus the root always has id equal to the number of nodes in the tree.
//
// Because this integer assignment cannot be known as we
// are building up the tree, we must make a second pass through the tree, which
// must mutate state. However, this re-id-ing pass is itself deterministic, so
// doing it a second time will always give the same result.
//
// leaves_ is a bitset indicating the set of leaves below. Similarly it needs to
// be calculated on a second pass, because we don't even know the size of the
// bitset as the tree is being built.
//
// Both of these features are prepared using the Polish method.
//
// In summary, call Polish after building your tree if you need to use internal
// node ids or leaf sets. Note that Tree construction calls Polish, if you are
// manually manipulating the topology make you do manipulations with that in
// mind.
//
// Equality is in terms of tree topologies. These mutable members don't matter.

#ifndef SRC_NODE_HPP_
#define SRC_NODE_HPP_

#include <limits.h>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "bitset.hpp"
#include "intpack.hpp"
#include "sugar.hpp"

class Node {
 public:
  typedef std::shared_ptr<Node> NodePtr;
  typedef std::vector<NodePtr> NodePtrVec;
  typedef std::shared_ptr<NodePtrVec> NodePtrVecPtr;
  typedef std::unordered_map<NodePtr, uint32_t> TopologyCounter;

  // This is the type of functions that are used in the PCSS recursion
  // functions. See `doc/pcss.svg` for a diagram of the PCSS traversal. In that
  // file, the first tree shows the terminology, and the subsequent trees show
  // the calls to f_root and f_internal.
  //
  // The signature is in 5 parts. The first 4 describe the position in the tree
  // and then the direction: the sister clade, the focal clade, child 0, and
  // child 1. False means down the tree structure and true means up. The 5th
  // part is the top of the virtual root clade, namely the clade containing the
  // virtual root (shown in gray in the diagram). Caution: in the case where the
  // virtual root clade is above the subsplit, the "virtual root clade" will be
  // the entire tree. There's nothing else we can do without rerooting the tree.
  // It's not too hard to exclude the undesired bits with a conditional tree
  // traversal. See IndexerRepresentationOfTopology for an example.
  using PCSSFun = std::function<void(const Node*, bool, const Node*, bool, const Node*,
                                     bool, const Node*, bool, const Node*)>;
  // The rooted version just uses: sister clade, the focal clade, child 0, and child 1.
  using RootedPCSSFun =
      std::function<void(const Node*, const Node*, const Node*, const Node*)>;

 public:
  explicit Node(uint32_t leaf_id, Bitset leaves);
  explicit Node(NodePtrVec children, size_t id, Bitset leaves);

  size_t Id() const { return id_; }
  uint64_t Tag() const { return tag_; }
  const Bitset& Leaves() const { return leaves_; }
  std::string TagString() const { return StringOfPackedInt(this->tag_); }
  uint32_t MaxLeafID() const { return MaxLeafIDOfTag(tag_); }
  uint32_t LeafCount() const { return LeafCountOfTag(tag_); }
  size_t Hash() const { return hash_; }
  bool IsLeaf() const { return children_.empty(); }
  const NodePtrVec& Children() const { return children_; }

  bool operator==(const Node& other) const;

  void PreOrder(std::function<void(const Node*)> f) const;
  // ConditionalPreOrder continues to recur as long as f returns true.
  void ConditionalPreOrder(std::function<bool(const Node*)> f) const;
  void PostOrder(std::function<void(const Node*)> f) const;
  void LevelOrder(std::function<void(const Node*)> f) const;
  // Apply the pre function before recurring down the tree, and then apply the
  // post function as we are recurring back up the tree.
  void PrePostOrder(std::function<void(const Node*)> pre,
                    std::function<void(const Node*)> post) const;

  // We take in two functions, f_root, and f_internal, each of which take three
  // edges.
  // We assume that f_root is symmetric in its last two arguments so that
  // f_root's signature actually looks like f_root(node0, {node1, node2}).
  // We apply f_root to the descendant edges like so: 012, 120, and 201. Because
  // f_root is symmetric in the last two arguments, we are going to get all of
  // the distinct calls of f.
  // At the internal nodes we cycle through triples of (node, sister, parent)
  // for f_internal.
  void TriplePreOrder(
      std::function<void(const Node*, const Node*, const Node*)> f_root,
      std::function<void(const Node*, const Node*, const Node*)> f_internal) const;
  // Iterate f through (node, sister, parent) for bifurcating trees using a
  // preorder traversal.
  void TriplePreOrderBifurcating(
      std::function<void(const Node*, const Node*, const Node*)> f) const;
  // As above, but getting indices rather than nodes themselves.
  void TripleIdPreOrderBifurcating(std::function<void(int, int, int)> f) const;

  // These two functions take functions accepting triples of (node_id,
  // child0_id, child1_id) and apply them according to various traversals.
  void BinaryIdPreOrder(const std::function<void(int, int, int)> f) const;
  void BinaryIdPostOrder(const std::function<void(int, int, int)> f) const;

  // See the typedef of PCSSFun and RootedPCSSFun to understand the argument type to
  // these functions.
  void PCSSPreOrder(PCSSFun f) const;
  void RootedPCSSPreOrder(RootedPCSSFun f) const;

  // This function prepares the id_ and leaves_ member variables as described at
  // the start of this document. It returns a map that maps the tags to their
  // indices. It's the verb, not the nationality.
  TagSizeMap Polish();

  // Return a vector such that the ith component describes the indices for nodes
  // above the current node.
  SizeVectorVector IdsAbove() const;

  std::string Newick(std::function<std::string(const Node*)> node_labeler,
                     const DoubleVectorOption& branch_lengths = std::nullopt) const;

  std::string Newick(const DoubleVectorOption& branch_lengths = std::nullopt,
                     const TagStringMapOption& node_labels = std::nullopt,
                     bool show_tags = false) const;

  // Construct a vector such that the ith entry is the id of the parent of the
  // node having id i. We assume that the indices are contiguous, and that the
  // root has the largest id.
  std::vector<size_t> ParentIdVector() const;

  NodePtr Deroot();

  // ** Static methods
  static inline uint32_t MaxLeafIDOfTag(uint64_t tag) { return UnpackFirstInt(tag); }
  static inline uint32_t LeafCountOfTag(uint64_t tag) { return UnpackSecondInt(tag); }
  static NodePtr Leaf(uint32_t id, Bitset leaves = Bitset(0));
  // Join builds a Node with the given descendants, or-ing the leaves_ of the
  // descendants.
  static NodePtr Join(NodePtrVec children, size_t id = SIZE_MAX);
  static NodePtr Join(NodePtr left, NodePtr right, size_t id = SIZE_MAX);
  // Build a tree given a vector of indices, such that each entry gives the
  // id of its parent. We assume that the indices are contiguous, and that
  // the root has the largest id.
  static NodePtr OfParentIdVector(std::vector<size_t> indices);

  //     topology           with internal node indices
  //     --------           --------------------------
  // 0: (0,1,(2,3))         (0,1,(2,3)4)5;
  // 1; (0,1,(2,3)) again   (0,1,(2,3)4)5;
  // 2: (0,2,(1,3))         (0,2,(1,3)4)5;
  // 3: (0,(1,(2,3)))       (0,(1,(2,3)4)5)6;
  static NodePtrVec ExampleTopologies();

  // Make a maximally-unbalanced "ladder" tree.
  static NodePtr Ladder(uint32_t leaf_count);

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
  // See beginning of file for notes about the id and the leaves.
  size_t id_;
  Bitset leaves_;
  // The tag_ is a pair of packed integers representing (1) the maximum leaf ID
  // of the leaves below this node, and (2) the number of leaves below the node.
  uint64_t tag_;
  size_t hash_;

  // Make copy constructors private to eliminate copying.
  Node(const Node&);
  Node& operator=(const Node&);

  // This is a private PostOrder that can change the Node.
  void MutablePostOrder(std::function<void(Node*)> f);

  std::string NewickAux(std::function<std::string(const Node*)> node_labeler,
                        const DoubleVectorOption& branch_lengths) const;

  // Make a leaf bitset by or-ing the leaf bitsets of the provided children.
  // Private just to avoid polluting the public interface.
  static Bitset LeavesOf(const NodePtrVec& children);
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

typedef std::unordered_map<uint64_t, Bitset> TagBitsetMap;

// Make a map from Tags to the bitset representing the leaves below the Tag.
// Just used for testing now.
TagBitsetMap TagLeafSetMapOf(Node::NodePtr topology) {
  TagBitsetMap map;
  auto leaf_count = topology->LeafCount();
  topology->PostOrder([&map, leaf_count](const Node* node) {
    Bitset bitset((size_t)leaf_count);
    if (node->IsLeaf()) {
      bitset.set(node->MaxLeafID());
    } else {
      // Take the union of the children below.
      for (const auto& child : node->Children()) {
        bitset |= map.at(child->Tag());
      }
    }
    SafeInsert(map, node->Tag(), std::move(bitset));
  });
  return map;
}

TEST_CASE("Node") {
  Node::NodePtrVec examples = Node::ExampleTopologies();
  Node::NodePtr t1 = examples[0];       // 0: (0,1,(2,3))
  Node::NodePtr t1_twin = examples[1];  // 1; (0,1,(2,3)) again
  Node::NodePtr t2 = examples[2];       // 2: (0,2,(1,3))
  Node::NodePtr t3 = examples[3];       // 3: (0,(1,(2,3)))
  // ((((0,1)7,2)8,(3,4)9)10,5,6)11;
  Node::NodePtr t4 = Node::OfParentIdVector({7, 7, 8, 9, 9, 11, 11, 8, 10, 10, 11});

  std::vector<std::string> triples;
  auto collect_triple = [&triples](const Node* node, const Node* sister,
                                   const Node* parent) {
    triples.push_back(std::to_string(node->Id()) + ", " + std::to_string(sister->Id()) +
                      ", " + std::to_string(parent->Id()));
  };
  t4->TriplePreOrder(collect_triple, collect_triple);
  std::vector<std::string> correct_triples(
      {"10, 5, 6", "8, 9, 10", "7, 2, 8", "0, 1, 7", "1, 0, 7", "2, 7, 8", "9, 8, 10",
       "3, 4, 9", "4, 3, 9", "5, 6, 10", "6, 10, 5"});
  CHECK_EQ(triples, correct_triples);

  // This is actually a non-trivial test (see note in Node constructor above),
  // which shows why we need bit rotation.
  CHECK_NE(t1->Hash(), t2->Hash());

  CHECK_EQ(t1, t1_twin);
  CHECK_NE(t1, t2);

  // Tree with trifurcation at the root.
  Node::NodePtr t1_alt = Node::OfParentIdVector({5, 5, 4, 4, 5});
  CHECK_EQ(t1, t1_alt);
  // Bifurcating tree.
  Node::NodePtr t3_alt = Node::OfParentIdVector({6, 5, 4, 4, 5, 6});
  CHECK_EQ(t3, t3_alt);

  for (const auto& topology : examples) {
    CHECK_EQ(topology, Node::OfParentIdVector(topology->ParentIdVector()));
    auto tag_leaf_set_map = TagLeafSetMapOf(topology);
    topology->PreOrder([&tag_leaf_set_map](const Node* node) {
      CHECK_EQ(node->Leaves(), tag_leaf_set_map.at(node->Tag()));
    });
  }

  // Check Deroot when we deroot on the right.
  CHECK_EQ(t1, t3->Deroot());
  // Check Deroot when we deroot on the left.
  CHECK_EQ(Node::OfParentIdVector({3, 3, 3}),
           // tree ((0,1)3,2)4
           Node::OfParentIdVector({3, 3, 4, 4})->Deroot());

  CHECK_EQ(Node::OfParentIdVector({4, 4, 5, 6, 5, 6}), Node::Ladder(4));
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_NODE_HPP_
