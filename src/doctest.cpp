#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "sbn.hpp"

TEST_CASE("Trying out Node") {
  auto t = Node::Join(Node::Join(Node::Leaf(0), Node::Leaf(1)), Node::Leaf(2));

  auto print_pos = [](Node* t) {
    std::cout << "I'm at " << t->TagString() << std::endl;
  };
  std::cout << t->ToNewick() << std::endl;
  std::cout << "PreOrder\n";
  t->PreOrder(print_pos);
  std::cout << "\nPostOrder\n";
  t->PostOrder(print_pos);

  // std::cout t->
}

