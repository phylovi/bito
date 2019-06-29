// TODO(erick)
// add branch lengths
// make everything const


#ifndef SRC_TREE_HPP_
#define SRC_TREE_HPP_

#include <limits.h>
#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "intpack.hpp"

class Node {
 public:
  typedef std::shared_ptr<Node> NodePtr;
  typedef std::vector<NodePtr> NodePtrVec;
  typedef std::shared_ptr<NodePtrVec> NodePtrVecPtr;
  typedef std::unordered_map<NodePtr, unsigned int> NodePtrCounter;
  typedef std::shared_ptr<NodePtrCounter> NodePtrCounterPtr;

 private:
  NodePtrVec children_;
  // The tag_ is a pair of packed integers representing (1) the maximum leaf ID
  // of the leaves below this node, and (2) the number of leaves below the node.
  uint64_t tag_;
  size_t hash_;

  // Make copy constructors private to eliminate copying.
  Node(const Node&);
  Node& operator=(const Node&);

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

  uint64_t Tag() { return tag_; };
  std::string TagString() { return StringOfPackedInt(this->tag_); }
  uint32_t MaxLeafID() const { return UnpackFirstInt(tag_); }
  uint32_t LeafCount() const { return UnpackSecondInt(tag_); }
  size_t Hash() const { return hash_; }
  bool IsLeaf() { return children_.empty(); }
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

  void PreOrder(std::function<void(Node*)> f) {
    f(this);
    for (auto child : children_) {
      child->PreOrder(f);
    }
  }

  // Specialized tree traversal for the standard representation of "unrooted"
  // trees, that have a trifurcation at the root.
  // Does not apply f to root node.
  // Node, Parent, Sister, parent_is_root
  void NPSPreOrderAux(std::function<void(Node*, Node*, Node*, bool)> f,
                      bool parent_is_root) {
    if (!IsLeaf()) {
      assert(children_.size() == 2);
      f(children_[0].get(), this, children_[1].get(), parent_is_root);
      children_[0]->NPSPreOrderAux(f, false);
      f(children_[1].get(), this, children_[0].get(), parent_is_root);
      children_[1]->NPSPreOrderAux(f, false);
    }
  }
  void NPSPreOrder(std::function<void(Node*, Node*, Node*, bool)> f) {
    assert(children_.size() == 3);
    for (auto child : children_) {
      child->NPSPreOrderAux(f, true);
    }
  }

  void PostOrder(std::function<void(Node*)> f) {
    for (auto child : children_) {
      child->PostOrder(f);
    }
    f(this);
  }

  void LevelOrder(std::function<void(Node*)> f) {
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

  std::string Newick() { return NewickAux() + ";"; }

  std::string NewickAux() {
    if (IsLeaf()) {
      return TagString();
    }
    std::string str = "(";
    for (auto iter = children_.begin(); iter != children_.end(); iter++) {
      if (iter != children_.begin()) {
        str.append(",");
      }
      str.append((*iter)->NewickAux());
    }
    str.append(")");
    str.append(TagString());
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
}


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

  // TODO add real test for NPSPreorder
  t3->NPSPreOrder(
      [](Node* node, Node* parent, Node* sister, bool parent_is_root) {
        std::cout << node->TagString() << ", " << parent->TagString() << ", "
                  << sister->TagString() << " " << parent_is_root << std::endl;
      });

  // This is actually a non-trivial test (see note in Node constructor above),
  // which shows why we need bit rotation.
  CHECK_NE(t1->Hash(), t2->Hash());

  CHECK_EQ(t1, t1_twin);
  CHECK_NE(t1, t2);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_HPP_
