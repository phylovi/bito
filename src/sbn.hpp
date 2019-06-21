// TODO
// document tag
// add branch lengths

// To discuss:
// look through abort-- how to handle?
// unique_ptr? Are we copying things in parsing?
// "NodeId" type rather than unsigned int?

#ifndef __SBN_HPP
#define __SBN_HPP

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

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
  Node(unsigned int leaf_id)
      : children_({}), max_leaf_id_(leaf_id), leaf_count_(1) {}
  // TODO here we are doing a copy of children?
  Node(NodePtrVec children) {
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
                  std::cout << "Tie observed between " << lhs->ToNewick()
                            << " and " << rhs->ToNewick() << std::endl;
                  abort();
                }
                return (difference < 0);
              });
    max_leaf_id_ = children_.back()->MaxLeafID();
    leaf_count_ = 0;
    for (auto child : children_) {
      leaf_count_ += child->LeafCount();
    }
    // Children are sorted by their max_leaf_id, so we can get the max by
    // looking at the last element.
  }
  ~Node() {
    // std::cout << "Destroying node " << TagString() << std::endl;
  }

  unsigned int MaxLeafID() const { return max_leaf_id_; }
  unsigned int LeafCount() const { return leaf_count_; }
  bool IsLeaf() { return children_.empty(); }

  std::string TagString() {
    return "<" + std::to_string(max_leaf_id_) + "," +
           std::to_string(leaf_count_) + ">";
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

  // void LevelOrder(std::function<void(Node*)> f) {
  //   for (auto child : children_) {
  //     child->PostOrder(f);
  //   }
  //   f(this);
  // }

  //     def _iter_descendants_levelorder(self, is_leaf_fn=None):
  //         """
  //         Iterate over all desdecendant nodes.
  //         """
  //         tovisit = deque([self])
  //         while len(tovisit)>0:
  //             node = tovisit.popleft()
  //             yield node
  //             if not is_leaf_fn or not is_leaf_fn(node):
  //                 tovisit.extend(node.children)
  //


  std::vector<unsigned int> MaxLeafTrace() {
    std::vector<unsigned int> trace;
    PreOrder([&trace](Node* node) { trace.push_back(node->MaxLeafID()); });
    return trace;
  }

  std::string ToNewick() {
    if (IsLeaf()) {
      return TagString();
    }
    std::string str = "(";
    for (auto iter = children_.begin(); iter != children_.end(); iter++) {
      if (iter != children_.begin()) {
        str.append(",");
      }
      str.append((*iter)->ToNewick());
    }
    str.append(")");
    str.append(TagString());
    return str;
  }

  // Class methods
  static NodePtr Leaf(int id) { return std::make_shared<Node>(id); }
  static NodePtr Join(NodePtrVec children) {
    return std::make_shared<Node>(children);
  };
  static NodePtr Join(NodePtr left, NodePtr right) {
    return Join(std::vector<NodePtr>({left, right}));
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#endif
#endif
