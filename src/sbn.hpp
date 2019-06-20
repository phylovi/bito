#ifndef __SBN_HPP
#define __SBN_HPP

#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

class MyClass {
 private:
  int m_i;

 public:
  void int_set(int i);

  int int_get();
};


// Class for a tree.
// Has nodes with unsigned integer ids.
// These ids have to increase as we go towards the root.
// TODO we don't check to make sure that the ids are different.
class Node {
  typedef std::shared_ptr<Node> NodePtr;
  typedef std::vector<NodePtr> NodePtrVec;

 private:
  NodePtrVec children_;
  // TODO: "NodeId" type rather than unsigned int?
  unsigned int id_;

  // Make copy constructors private to eliminate copying.
  Node(const Node&);
  Node& operator=(const Node&);


 public:
  Node(unsigned int id) : children_({}), id_(id) {}
  Node(NodePtrVec children, unsigned int id) : children_(children), id_(id) {
    // Nodes must have a larger index than their children.
    assert(MaxChildIdx(children) < id);
  }
  Node(NodePtr left, NodePtr right, unsigned int id)
      : Node({left, right}, id) {}
  ~Node() { std::cout << "Destroying node " << id_ << std::endl; }

  unsigned int GetId() { return id_; }
  bool IsLeaf() { return children_.empty(); }

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

  unsigned int LeafCount() {
    unsigned int count = 0;
    PreOrder([&count](Node* node) { count += node->IsLeaf(); });
    return count;
  }

  // Class methods
  static NodePtr Leaf(int id) { return std::make_shared<Node>(id); }
  static NodePtr Join(NodePtr left, NodePtr right, int id) {
    return std::make_shared<Node>(left, right, id);
  };
  static unsigned int MaxChildIdx(NodePtrVec children) {
    assert(~children.empty());
    // 0 is the smallest value for an unsigned integer.
    unsigned int result = 0;
    for (auto child : children) {
      result = std::max(child->GetId(), result);
    }
    return result;
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Trying out Node") {
  auto t =
      Node::Join(Node::Join(Node::Leaf(0), Node::Leaf(1), 3), Node::Leaf(2), 4);

  REQUIRE(t->LeafCount() == 3);

  auto print_pos = [](Node* t) {
    std::cout << "I'm at " << t->GetId() << std::endl;
  };
  t->PreOrder(print_pos);
  t->PostOrder(print_pos);
}
#endif
#endif
