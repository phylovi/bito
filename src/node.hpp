// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_NODE_HPP_
#define SRC_NODE_HPP_

#include <limits.h>
#include <algorithm>
#include <cassert>
#include <deque>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
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

  // This is the type of functions that are used in the PCSS recursion
  // functions. The signature is in 4 parts, each of which describes the
  // position in the tree and then the direction: false means down the tree
  // structure and true means up. The 4 parts are the uncut parent, the cut
  // parent, child 0, and child 1.
  typedef std::function<void(const Node*, bool, const Node*, bool, const Node*,
                             bool, const Node*, bool)>
      PCSSFun;

 public:
  explicit Node(unsigned int leaf_id)
      : children_({}), tag_(PackInts(leaf_id, 1)), hash_(SOHash(leaf_id)) {}
  explicit Node(NodePtrVec children) {
    children_ = children;
    if (children_.empty()) {
      // This constructor is for internal nodes, so we can't allow children to
      // be empty.
      abort();
    }
    // Order the children by their max leaf ids.
    std::sort(children_.begin(), children_.end(),
              [](const auto& lhs, const auto& rhs) {
                int difference = lhs->MaxLeafID() - rhs->MaxLeafID();
                // Children should have non-overlapping leaf sets, so there
                // should not be ties.
                if (difference == 0) {
                  std::cout << "Tie observed between " << lhs->Newick()
                            << " and " << rhs->Newick() << std::endl;
                  std::cout << "Do you have a taxon name repeated?\n";
                  abort();
                }
                return (difference < 0);
              });
    // Children are sorted by their max_leaf_id, so we can get the max by
    // looking at the last element.
    uint32_t max_leaf_id = children_.back()->MaxLeafID();
    uint32_t leaf_count = 0;
    hash_ = 0;
    for (auto child : children_) {
      leaf_count += child->LeafCount();
      hash_ ^= child->Hash();
    }
    tag_ = PackInts(max_leaf_id, leaf_count);
    // Bit rotation is necessary because if we only XOR then we can get
    // collisions when identical tips are in different
    // ordered subtrees (an example is in below doctest).
    hash_ = SORotate(hash_, 1);
  }

  uint64_t Tag() const { return tag_; }
  std::string TagString() const { return StringOfPackedInt(this->tag_); }
  uint32_t MaxLeafID() const { return UnpackFirstInt(tag_); }
  uint32_t LeafCount() const { return UnpackSecondInt(tag_); }
  size_t Hash() const { return hash_; }
  bool IsLeaf() const { return children_.empty(); }
  NodePtrVec Children() const { return children_; }

  bool operator==(const Node& other) {
    if (this->Hash() != other.Hash()) {
      return false;
    }
    size_t child_count = this->Children().size();
    if (child_count != other.Children().size()) {
      return false;
    }
    for (size_t i = 0; i < child_count; i++) {
      if (!(*children_[i] == *other.Children()[i])) {
        return false;
      }
    }
    return true;
  }

  void PreOrder(std::function<void(const Node*)> f) {
    f(this);
    for (auto child : children_) {
      child->PreOrder(f);
    }
  }

  // Iterate f through (parent, sister, node) for internal nodes using a
  // preorder traversal.
  void TriplePreOrderInternal(
      std::function<void(const Node*, const Node*, const Node*)> f) {
    if (!IsLeaf()) {
      assert(children_.size() == 2);
      f(this, children_[1].get(), children_[0].get());
      children_[0]->TriplePreOrderInternal(f);
      f(this, children_[0].get(), children_[1].get());
      children_[1]->TriplePreOrderInternal(f);
    }
  }

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
      std::function<void(const Node*, const Node*, const Node*)> f_internal) {
    assert(children_.size() == 3);
    f_root(children_[0].get(), children_[1].get(), children_[2].get());
    f_root(children_[1].get(), children_[2].get(), children_[0].get());
    f_root(children_[2].get(), children_[0].get(), children_[1].get());
    for (auto child : children_) {
      child->TriplePreOrderInternal(f_internal);
    }
  }

  // See the typedef of PCSSFun to understand the argument type to this
  // function.
  void PCSSPreOrder(PCSSFun f) {
    this->TriplePreOrder(
        // f_root
        [&f](const Node* node0, const Node* node1, const Node* node2) {
          // Virtual root on node2's edge, with subsplit pointing up.
          f(node2, false, node2, true, node0, false, node1, false);
          if (!node2->IsLeaf()) {
            assert(node2->Children().size() == 2);
            auto child0 = node2->Children()[0].get();
            auto child1 = node2->Children()[1].get();
            // Virtual root in node1.
            f(node0, false, node2, false, child0, false, child1, false);
            // Virtual root in node0.
            f(node1, false, node2, false, child0, false, child1, false);
            // Virtual root on node2's edge, with subsplit pointing down.
            f(node2, true, node2, false, child0, false, child1, false);
            // Virtual root in child0.
            f(child1, false, node2, true, node0, false, node1, false);
            // Virtual root in child1.
            f(child0, false, node2, true, node0, false, node1, false);
          }
        },
        // f_internal
        [&f](const Node* parent, const Node* sister, const Node* node) {
          // Virtual root on node's edge, with subsplit pointing up.
          f(node, false, node, true, parent, true, sister, false);
          if (!node->IsLeaf()) {
            assert(node->Children().size() == 2);
            auto child0 = node->Children()[0].get();
            auto child1 = node->Children()[1].get();
            // Virtual root up the tree.
            f(sister, false, node, false, child0, false, child1, false);
            // Virtual root in sister.
            f(parent, true, node, false, child0, false, child1, false);
            // Virtual root on node's edge, with subsplit pointing down.
            f(node, true, node, false, child0, false, child1, false);
            // Virtual root in child0.
            f(child1, false, node, true, sister, false, parent, true);
            // Virtual root in child1.
            f(child0, false, node, true, sister, false, parent, true);
          }
        });
  }

  void PostOrder(std::function<void(const Node*)> f) {
    for (auto child : children_) {
      child->PostOrder(f);
    }
    f(this);
  }

  int PostOrder(std::function<int(const Node*, const std::vector<int>&)> f) {
    std::vector<int> v;
    for (auto child : children_) {
      v.push_back(child->PostOrder(f));
    }
    return f(this, v);
  }

  void LevelOrder(std::function<void(const Node*)> f) {
    std::deque<Node*> to_visit = {this};
    while (to_visit.size()) {
      auto n = to_visit.front();
      f(n);
      to_visit.pop_front();

      for (auto child : n->children_) {
        to_visit.push_back(child.get());
      }
    }
  }

  std::string Newick(
      const TagDoubleMapOption& branch_lengths = std::experimental::nullopt,
      const TagStringMapOption& node_labels = std::experimental::nullopt) {
    return NewickAux(branch_lengths, node_labels) + ";";
  }

  std::string NewickAux(const TagDoubleMapOption& branch_lengths,
                        const TagStringMapOption& node_labels) {
    std::string str;
    if (IsLeaf()) {
      if (node_labels) {
        str.assign((*node_labels).at(Tag()));
      } else {
        str.assign(TagString());
      }
    } else {
      str.assign("(");
      for (auto iter = children_.begin(); iter != children_.end(); iter++) {
        if (iter != children_.begin()) {
          str.append(",");
        }
        str.append((*iter)->NewickAux(branch_lengths, node_labels));
      }
      str.append(")");
      if (!node_labels) {
        // If node_labels are not included then we figure that the discrete
        // structure of the tree is of interest. Thus we write out the tags as
        // internal nodes of the tree.
        str.append(TagString());
      }
    }
    if (branch_lengths) {
      auto search = (*branch_lengths).find(Tag());
      if (search != (*branch_lengths).end()) {
        // ostringstream is the way to get scientific notation using the STL.
        std::ostringstream str_stream;
        str_stream << search->second;
        str.append(":" + str_stream.str());
      }
    }
    return str;
  }

  // Class methods
  static NodePtr Leaf(int id) { return std::make_shared<Node>(id); }
  static NodePtr Join(NodePtrVec children) {
    return std::make_shared<Node>(children);
  }
  static NodePtr Join(NodePtr left, NodePtr right) {
    return Join(std::vector<NodePtr>({left, right}));
  }

  // A "cryptographic" hash function from Stack Overflow (the std::hash function
  // appears to leave unsigned ints as they are, which doesn't work for our
  // application).
  // https://stackoverflow.com/a/12996028/467327
  static inline unsigned int SOHash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
  }

  // Bit rotation from Stack Overflow.
  // c is the amount by which we rotate.
  // https://stackoverflow.com/a/776523/467327
  static inline size_t SORotate(size_t n, unsigned int c) {
    const unsigned int mask =
        (CHAR_BIT * sizeof(n) - 1);  // assumes width is a power of 2.
    // assert ( (c<=mask) &&"rotate by type width or more");
    c &= mask;
    return (n << c) | (n >> ((-c) & mask));
  }

 private:
  NodePtrVec children_;
  // The tag_ is a pair of packed integers representing (1) the maximum leaf ID
  // of the leaves below this node, and (2) the number of leaves below the node.
  uint64_t tag_;
  size_t hash_;

  // Make copy constructors private to eliminate copying.
  Node(const Node&);
  Node& operator=(const Node&);
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
  auto t1 = Node::Join(
      std::vector<Node::NodePtr>({Node::Leaf(0), Node::Leaf(1),
                                  Node::Join(Node::Leaf(2), Node::Leaf(3))}));
  auto t1_twin = Node::Join(
      std::vector<Node::NodePtr>({Node::Leaf(1), Node::Leaf(0),
                                  Node::Join(Node::Leaf(3), Node::Leaf(2))}));
  auto t2 = Node::Join(
      std::vector<Node::NodePtr>({Node::Leaf(0), Node::Leaf(2),
                                  Node::Join(Node::Leaf(1), Node::Leaf(3))}));
  auto t3 = Node::Join(std::vector<Node::NodePtr>(
      {Node::Leaf(0), Node::Leaf(1),
       Node::Join(Node::Leaf(2), Node::Join(Node::Leaf(3), Node::Leaf(4)))}));

  // TODO(ematsen) add real test for TriplePreorder
  std::cout << "TriplePreOrder" << std::endl;
  std::cout << t3->Newick() << std::endl;
  auto print_triple = [](const Node* parent, const Node* sister,
                         const Node* node) {
    std::cout << parent->TagString() << ", " << sister->TagString() << ", "
              << node->TagString() << " " << std::endl;
  };
  t3->TriplePreOrder(print_triple, print_triple);

  // This is actually a non-trivial test (see note in Node constructor above),
  // which shows why we need bit rotation.
  CHECK_NE(t1->Hash(), t2->Hash());

  CHECK_EQ(t1, t1_twin);
  CHECK_NE(t1, t2);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_NODE_HPP_
