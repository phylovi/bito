#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "sbn.hpp"

TEST_CASE("Trying out Node") {
  auto t = Node::Join(Node::Join(Node::Leaf(0), Node::Leaf(1)), Node::Leaf(2));

  auto print_pos = [](Node* t) {
    std::cout << "I'm at " << t->TagString() << std::endl;
  };
  std::cout << t->Newick() << std::endl;
  std::cout << "PreOrder\n";
  t->PreOrder(print_pos);
  std::cout << "\nPostOrder\n";
  t->PostOrder(print_pos);


  std::vector<std::string> trace;
  t->PreOrder([&trace](Node* node) { trace.push_back(node->TagString()); });

  for (auto x : trace) {
    std::cout << x << " ";
  }


  // preorder:
  REQUIRE(std::vector<std::string>({"9_10", "4_5", "3_4", "1_2", "0_1", "1_1",
                                    "3_2", "2_1", "3_1", "4_1", "9_5", "7_3",
                                    "5_1", "7_2", "6_1", "7_1", "9_2", "8_1",
                                    "9_1"}) == trace);


  // postorder:
  // REQUIRE(
  //     std::vector<std::string>({"0_1", "1_1", "1_2", "2_1", "3_1", "3_2",
  //     "3_4",
  //                               "4_1", "4_5", "5_1", "6_1", "7_1", "7_2",
  //                               "7_3",
  //                               "8_1", "9_1", "9_2", "9_5", "9_10"}) ==
  //                               trace);

  // levelorder:
  // REQUIRE(std::vector<std::string>({"9_10","4_5","9_5","3_4","4_1","7_3","9_2","1_2","3_2","5_1","7_2","8_1","9_1","0_1","1_1","2_1","3_1","6_1","7_1"})
  // == trace);
}

