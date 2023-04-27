#pragma once

#include "../src/node.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED

typedef std::unordered_map<uint64_t, Bitset> TagBitsetMap;

// Make a map from Tags to the bitset representing the leaves below the Tag.
// Just used for testing now.
TagBitsetMap TagLeafSetMapOf(Node::NodePtr topology) {
  TagBitsetMap map;
  auto leaf_count = topology->LeafCount();
  topology->Postorder([&map, leaf_count](const Node* node) {
    Bitset bitset(static_cast<size_t>(leaf_count));
    if (node->IsLeaf()) {
      bitset.set(node->MaxLeafID());
    } else {
      // Take the union of the children below.
      for (const auto& child : node->Children()) {
        bitset |= map.at(child->Tag());
      }
    }
    SafeInsert(map, node->Tag(), std::move(bitset));
  });
  return map;
}

TEST_CASE("Node") {
  Node::NodePtrVec examples = Node::ExampleTopologies();
  Node::NodePtr t1 = examples[0];       // 0: (0,1,(2,3))
  Node::NodePtr t1_twin = examples[1];  // 1; (0,1,(2,3)) again
  Node::NodePtr t2 = examples[2];       // 2: (0,2,(1,3))
  Node::NodePtr t3 = examples[3];       // 3: (0,(1,(2,3)))
  // ((((0,1)7,2)8,(3,4)9)10,5,6)11;
  Node::NodePtr tbig = Node::OfParentIdVector({7, 7, 8, 9, 9, 11, 11, 8, 10, 10, 11});

  std::vector<std::string> triples;
  auto collect_triple = [&triples](const Node* node, const Node* sister,
                                   const Node* parent) {
    triples.push_back(std::to_string(node->Id()) + ", " + std::to_string(sister->Id()) +
                      ", " + std::to_string(parent->Id()));
  };
  tbig->TriplePreorder(collect_triple, collect_triple);
  std::vector<std::string> correct_triples(
      {"10, 5, 6", "8, 9, 10", "7, 2, 8", "0, 1, 7", "1, 0, 7", "2, 7, 8", "9, 8, 10",
       "3, 4, 9", "4, 3, 9", "5, 6, 10", "6, 10, 5"});
  CHECK_EQ(triples, correct_triples);

  // This is actually a non-trivial test (see note in Node constructor above),
  // which shows why we need bit rotation.
  CHECK_NE(t1->Hash(), t2->Hash());

  CHECK_EQ(t1, t1_twin);
  CHECK_NE(t1, t2);

  // Tree with trifurcation at the root.
  Node::NodePtr t1_alt = Node::OfParentIdVector({5, 5, 4, 4, 5});
  CHECK_EQ(t1, t1_alt);
  // Bifurcating tree.
  Node::NodePtr t3_alt = Node::OfParentIdVector({6, 5, 4, 4, 5, 6});
  CHECK_EQ(t3, t3_alt);

  for (const auto& topology : examples) {
    CHECK_EQ(topology, Node::OfParentIdVector(topology->ParentIdVector()));
    CHECK_EQ(topology, topology->DeepCopy());
    auto tag_leaf_set_map = TagLeafSetMapOf(topology);
    topology->Preorder([&tag_leaf_set_map](const Node* node) {
      CHECK_EQ(node->Leaves(), tag_leaf_set_map.at(node->Tag()));
    });
  }

  // Check Deroot when we deroot on the right.
  CHECK_EQ(t1, t3->Deroot());
  // Check Deroot when we deroot on the left.
  CHECK_EQ(Node::OfParentIdVector({3, 3, 3}),
           // tree ((0,1)3,2)4
           Node::OfParentIdVector({3, 3, 4, 4})->Deroot());

  CHECK_EQ(Node::OfParentIdVector({4, 4, 5, 6, 5, 6}), Node::Ladder(4));

  SizeVector correct_sisters({5, 4, 3, 2});
  SizeVector sisters;
  t3->RootedSisterAndLeafTraversal([&sisters](const Node* sister, const Node* leaf) {
    sisters.push_back(sister->Id());
  });
  CHECK_EQ(correct_sisters, sisters);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
