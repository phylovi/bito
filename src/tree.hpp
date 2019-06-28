// TODO(erick)
// add branch lengths
// make everything const


#ifndef SRC_TREE_HPP_
#define SRC_TREE_HPP_

#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "intpack.hpp"

class Node {
 public:
  typedef std::shared_ptr<Node> NodePtr;
  typedef std::vector<NodePtr> NodePtrVec;
  typedef std::shared_ptr<NodePtrVec> NodePtrVecPtr;

 private:
  NodePtrVec children_;
  // The tag_ is a pair of packed integers representing (1) the maximum leaf ID
  // of the leaves below this node, and (2) the number of leaves below the node.
  uint64_t tag_;

  // Make copy constructors private to eliminate copying.
  Node(const Node&);
  Node& operator=(const Node&);

 public:
  explicit Node(unsigned int leaf_id)
      : children_({}), tag_(PackInts(leaf_id, 1)) {}
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
    for (auto child : children_) {
      leaf_count += child->LeafCount();
    }
    tag_ = PackInts(max_leaf_id, leaf_count);
  }

  uint64_t Tag() { return tag_; };
  uint32_t MaxLeafID() const { return UnpackFirstInt(tag_); }
  uint32_t LeafCount() const { return UnpackSecondInt(tag_); }
  bool IsLeaf() { return children_.empty(); }
  NodePtrVec Children() const { return children_; }

  std::string TagString() { return StringOfPackedInt(this->tag_); }

  void PreOrder(std::function<void(Node*)> f) {
    f(this);
    for (auto child : children_) {
      child->PreOrder(f);
    }
  }

  void NPSPreOrderAux(std::function<void(Node*, Node*, Node*)> f) {
    if (!IsLeaf()) {
      assert(children_.size() == 2);
      f(children_[0].get(), this, children_[1].get());
      children_[0]->NPSPreOrderAux(f);
      f(children_[1].get(), this, children_[0].get());
      children_[1]->NPSPreOrderAux(f);
    }
  }
  void NPSPreOrder(std::function<void(Node*, Node*, Node*)> f) {
    for (auto child : children_) {
      child->NPSPreOrderAux(f);
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

  std::vector<unsigned int> MaxLeafTrace() {
    std::vector<unsigned int> trace;
    PreOrder([&trace](Node* node) { trace.push_back(node->MaxLeafID()); });
    return trace;
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
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Node header") {
  auto t = Node::Join(
      std::vector<Node::NodePtr>({Node::Leaf(0), Node::Leaf(1),
                                  Node::Join(Node::Leaf(2), Node::Leaf(3))}));

  // TODO add real test for NPSPreorder
  t->NPSPreOrder([](Node* node, Node* parent, Node* sister) {
    std::cout << node->TagString() << ", " << parent->TagString() << ", "
              << sister->TagString() << std::endl;
  });
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_HPP_
