#ifndef __SBN_HPP
#define __SBN_HPP

#include "doctest.h"

#include <memory>
#include <vector>

class MyClass {
 private:
  int m_i;

 public:
  void int_set(int i);

  int int_get();
};


class Node {

 typedef std::shared_ptr<Node> NodePtr;
 typedef std::vector<NodePtr> NodePtrVec;

 private:
  NodePtrVec children_;
  int id_;


 public:
  Node(NodePtrVec children, int id) : children_(children), id_(id) {}

  Node(NodePtr left, NodePtr right, int id) {
    children_ = {left, right};
    id_ = id;
  }

  int get_id() { return id_; }
};


TEST_CASE("testing the factorial function") {
    CHECK(1 == 1);
}

#endif
