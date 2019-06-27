// TODO(erick)
// document tag
// add branch lengths

// To discuss:
// look through abort-- how to handle? also look for cassert
// unique_ptr? Are we copying things in parsing?
// "NodeId" type rather than unsigned int?

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
  unsigned int max_leaf_id_;
  unsigned int leaf_count_;

  // Make copy constructors private to eliminate copying.
  Node(const Node&);
  Node& operator=(const Node&);

 public:
  explicit Node(unsigned int leaf_id)
      : children_({}), max_leaf_id_(leaf_id), leaf_count_(1) {}
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
    max_leaf_id_ = children_.back()->MaxLeafID();
    leaf_count_ = 0;
    for (auto child : children_) {
      leaf_count_ += child->LeafCount();
    }
  }

  unsigned int MaxLeafID() const { return max_leaf_id_; }
  unsigned int LeafCount() const { return leaf_count_; }
  bool IsLeaf() { return children_.empty(); }

  std::string TagString() {
    return std::to_string(max_leaf_id_) + "_" + std::to_string(leaf_count_);
  }

  void PreOrder(std::function<void(Node*)> f) {
    f(this);
    for (auto child : children_) {
      child->PreOrder(f);
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

#endif  // SRC_TREE_HPP_
