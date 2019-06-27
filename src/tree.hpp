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

  std::string Newick() {
    std::function<std::string(NodePtr)> Aux;
    Aux = [&Aux](NodePtr n) {
      if (n->IsLeaf()) {
        return n->TagString();
      }
      std::string str = "(";
      for (auto iter = n->Children().begin(); iter != n->Children().end();
           ++iter) {
        if (iter != n->Children().begin()) {
          str.append(",");
        }
        str.append(Aux(*iter));
      }
      str.append(")");
      str.append(n->TagString());
      return str;
    };
    std::shared_ptr<Node> this_shared(this);
    return Aux(this_shared) + ";";
  }

  // Class methods
  static NodePtr
  Leaf(int id) {
    return std::make_shared<Node>(id);
  }
  static NodePtr Join(NodePtrVec children) {
    return std::make_shared<Node>(children);
  }
  static NodePtr Join(NodePtr left, NodePtr right) {
    return Join(std::vector<NodePtr>({left, right}));
  }
};

#endif  // SRC_TREE_HPP_
