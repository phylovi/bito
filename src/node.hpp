// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// The Node class is how we express tree topologies.
//
// Nodes are immutable after construction except for the id_. The id_ is
// provided for applications where it is useful to have the edges numbered with
// a contiguous set of integers, where the leaf edges have the smallest such
// integers. Because this integer assignment cannot be known as we are building
// up the tree, we must make a second pass through the tree, which must mutate
// state. However, this re-id-ing pass is itself deterministic, so doing it a
// second time will always give the same result.
//
// In summary, call Reid after building your tree if you need to use internal
// node ids. Note that Tree construction calls Reid, if you are manually
// manipulating the topology make you do manipulations with that in mind.
//
// Equality is in terms of tree topologies. Indices don't matter.

#ifndef SRC_NODE_HPP_
#define SRC_NODE_HPP_

#include <limits.h>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
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
  typedef std::function<void(const Node*, bool, const Node*, bool, const Node*,
                             bool, const Node*, bool, const Node*)>
      PCSSFun;

 public:
  explicit Node(uint32_t leaf_id);
  explicit Node(NodePtrVec children, size_t id);

  size_t Index() const { return id_; }
  uint64_t Tag() const { return tag_; }
  std::string TagString() const { return StringOfPackedInt(this->tag_); }
  uint32_t MaxLeafID() const { return MaxLeafIDOfTag(tag_); }
  uint32_t LeafCount() const { return LeafCountOfTag(tag_); }
  size_t Hash() const { return hash_; }
  bool IsLeaf() const { return children_.empty(); }
  NodePtrVec Children() const { return children_; }

  bool operator==(const Node& other);

  void PreOrder(std::function<void(const Node*)> f) const;
  // ConditionalPreOrder continues to recur as long as f returns true.
  void ConditionalPreOrder(std::function<bool(const Node*)> f) const;
  void PostOrder(std::function<void(const Node*)> f) const;
  void LevelOrder(std::function<void(const Node*)> f) const;
  // Apply the pre function before recurring down the tree, and then apply the
  // post function as we are recurring back up the tree.
  void PrePostOrder(std::function<void(const Node*)> pre,
                    std::function<void(const Node*)> post) const;

  // Iterate f through (parent, sister, node) for bifurcating trees using a
  // preorder traversal.
  void TriplePreOrderBifurcating(
      std::function<void(const Node*, const Node*, const Node*)> f) const;
  // As above, but getting indices rather than nodes themselves.
  void TripleIndexPreOrderBifurcating(
      std::function<void(int, int, int)> f) const;

  // These two functions take functions accepting triples of (node_id,
  // child0_id, child1_id) and apply them according to various traversals.
  void BinaryIndexPreOrder(const std::function<void(int, int, int)> f) const;
  void BinaryIndexPostOrder(const std::function<void(int, int, int)> f) const;

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
  // This one is just the part of the above that's run for the internal nodes.
  void TriplePreOrderInternal(
      std::function<void(const Node*, const Node*, const Node*)> f) const;

  // See the typedef of PCSSFun to understand the argument type to this
  // function.
  void PCSSPreOrder(PCSSFun f) const;

  // This function assigns indices to the nodes of the topology: the leaves get
  // their indices (which are contiguously numbered from 0 through the leaf
  // count minus 1) and the rest get ordered according to a postorder traversal.
  // Thus the root always has id equal to the number of nodes in the tree.
  // It returns a map that maps the tags to their indices.
  TagSizeMap Reid();

  // Return a vector such that the ith component describes the indices for nodes
  // above the current node.
  SizeVectorVector IndicesAbove();

  std::string Newick(std::function<std::string(const Node*)> node_labeler,
                     const DoubleVectorOption& branch_lengths =
                         std::experimental::nullopt) const;

  std::string Newick(
      const DoubleVectorOption& branch_lengths = std::experimental::nullopt,
      const TagStringMapOption& node_labels = std::experimental::nullopt,
      bool show_tags = false) const;

  // Construct a vector such that the ith entry is the id of the parent of the
  // node having id i. We assume that the indices are contiguous, and that the
  // root has the largest id.
  std::vector<size_t> ParentIndexVector();

  NodePtr Deroot();

  // ** Class methods
  static inline uint32_t MaxLeafIDOfTag(uint64_t tag) {
    return UnpackFirstInt(tag);
  }
  static inline uint32_t LeafCountOfTag(uint64_t tag) {
    return UnpackSecondInt(tag);
  }
  static NodePtr Leaf(uint32_t id);
  static NodePtr Join(NodePtrVec children, size_t id = SIZE_MAX);
  static NodePtr Join(NodePtr left, NodePtr right, size_t id = SIZE_MAX);
  // Build a tree given a vector of indices, such that each entry gives the
  // id of its parent. We assume that the indices are contiguous, and that
  // the root has the largest id.
  static NodePtr OfParentIndexVector(std::vector<size_t> indices);

  //     topology           with internal node indices
  //     --------           --------------------------
  // 0: (0,1,(2,3))         (0,1,(2,3)4)5;
  // 1; (0,1,(2,3)) again   (0,1,(2,3)4)5;
  // 2: (0,2,(1,3))         (0,2,(1,3)4)5;
  // 3: (0,(1,(2,3)))       (0,(1,(2,3)4)5)6;
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
  // See beginning of file for notes about the id.
  size_t id_;
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
TEST_CASE("Node") {
  Node::NodePtrVec examples = Node::ExampleTopologies();
  Node::NodePtr t1 = examples[0];       // 0: (0,1,(2,3))
  Node::NodePtr t1_twin = examples[1];  // 1; (0,1,(2,3)) again
  Node::NodePtr t2 = examples[2];       // 2: (0,2,(1,3))
  Node::NodePtr t3 = examples[3];       // 3: (0,(1,(2,3)))
  // TODO(ematsen) add real test for TriplePreorder
  std::vector<std::string> triples;
  auto collect_triple = [&triples](const Node* parent, const Node* sister,
                                   const Node* node) {
    triples.push_back(parent->TagString() + ", " + sister->TagString() + ", " +
                      node->TagString());
  };
  t2->TriplePreOrder(collect_triple, collect_triple);
  // t2 with tags: (0_1,2_1,(1_1,3_1)3_2)
  std::vector<std::string> correct_triples({"0_1, 2_1, 3_2", "2_1, 3_2, 0_1",
                                            "3_2, 0_1, 2_1", "3_2, 3_1, 1_1",
                                            "3_2, 1_1, 3_1"});
  CHECK_EQ(triples, correct_triples);

  // This is actually a non-trivial test (see note in Node constructor above),
  // which shows why we need bit rotation.
  CHECK_NE(t1->Hash(), t2->Hash());

  CHECK_EQ(t1, t1_twin);
  CHECK_NE(t1, t2);

  // Tree with trifurcation at the root.
  Node::NodePtr t1_alt = Node::OfParentIndexVector({5, 5, 4, 4, 5});
  CHECK_EQ(t1, t1_alt);
  // Bifurcating tree.
  Node::NodePtr t3_alt = Node::OfParentIndexVector({6, 5, 4, 4, 5, 6});
  CHECK_EQ(t3, t3_alt);

  for (const auto& topology : examples) {
    CHECK_EQ(topology,
             Node::OfParentIndexVector(topology->ParentIndexVector()));
  }

  // Check Deroot when we deroot on the right.
  CHECK_EQ(t1, t3->Deroot());
  // Check Deroot when we deroot on the left.
  CHECK_EQ(Node::OfParentIndexVector({3, 3, 3}),
           // tree ((0,1)3,2)4
           Node::OfParentIndexVector({3, 3, 4, 4})->Deroot());
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_NODE_HPP_
