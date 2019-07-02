#ifndef SRC_BUILD_HPP_
#define SRC_BUILD_HPP_

#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "bitset.hpp"
#include "default_dict.hpp"
#include "driver.hpp"
#include "tree.hpp"

typedef std::unordered_map<uint64_t, Bitset> TagBitsetMap;
typedef std::unordered_set<Bitset> BitsetSet;
typedef std::unordered_map<Bitset, int> BitsetIndexer;
// TODO rename?
typedef DefaultDict<Bitset, uint32_t> BitsetUInt32Map;

// Using insert and at avoids needing to make a default constructor.
// https://stackoverflow.com/questions/17172080/insert-vs-emplace-vs-operator-in-c-map

// TODO do we want optional minor here?
TagBitsetMap TagBitsetMapOf(Node::NodePtr t) {
  TagBitsetMap m;
  auto leaf_count = t->LeafCount();
  t->PostOrder([&m, leaf_count](Node* n) {
    Bitset x((size_t)leaf_count);
    if (n->IsLeaf()) {
      x.set(n->MaxLeafID());
    } else {
      // Take the union of the children below.
      for (auto child : n->Children()) {
        x |= m.at(child->Tag());
      }
    }
    assert(m.insert(std::make_pair(n->Tag(), std::move(x))).second);
  });
  return m;
}

void PrintTagBitsetMap(TagBitsetMap m) {
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    std::cout << StringOfPackedInt(iter->first) << " "
              << iter->second.ToString() << std::endl;
  }
}

BitsetUInt32Map RootsplitCounterOf(Node::NodePtrCounterPtr trees) {
  BitsetUInt32Map rootsplit_counter(0);
  for (auto iter = trees->begin(); iter != trees->end(); ++iter) {
    auto tree = iter->first;
    auto count = iter->second;
    auto tag_to_bitset = TagBitsetMapOf(tree);
    auto Aux = [&rootsplit_counter, &tag_to_bitset, &count](Node* n) {
      rootsplit_counter.increment(tag_to_bitset.at(n->Tag()), count);
    };
    for (auto child : tree->Children()) {
      child->PreOrder(Aux);
    }
  }
  return rootsplit_counter;
}

// BitsetUInt32Map SupportsOf(Node::NodePtrCounterPtr trees) {
//   BitsetUInt32Map rootsplit_counter(0);
//   for (auto iter = trees->begin(); iter != trees->end(); ++iter) {
//     auto tree = iter->first;
//     auto count = iter->second;
//     auto tag_to_bitset = TagBitsetMapOf(tree);
//     tree->TriplePreOrder(
//         [&rootsplit_counter, &tag_to_bitset, &count](Node* n, Node*, Node*) {
//           rootsplit_counter.increment(tag_to_bitset.at(n->Tag()), count);
//         });
//   }
//   return rootsplit_counter;
// }

BitsetUInt32Map SupportsOf(Node::NodePtrCounterPtr trees) {
  BitsetUInt32Map subsplit_support(0);
  for (auto iter = trees->begin(); iter != trees->end(); ++iter) {
    auto tree = iter->first;
    auto count = iter->second;
    auto tag_to_bitset = TagBitsetMapOf(tree);
    auto leaf_count = tree->LeafCount();
    // TODO make a more informative error message when people don't put in a
    // bifurcating tree with a trifurcation at the root.
    tree->PCSSPreOrder([&subsplit_support, &tag_to_bitset, &count, &leaf_count](
        Node* parent_uncut_node, bool parent_uncut_direction,
        Node* parent_split_node, bool parent_cut_direction,  //
        Node* child0_node, bool child0_direction,            //
        Node* child1_node, bool child1_direction) {
      Bitset bitset(3 * leaf_count, false);
      bitset.CopyFrom(tag_to_bitset.at(parent_uncut_node->Tag()), 0,
                      parent_uncut_direction);
      bitset.CopyFrom(tag_to_bitset.at(parent_split_node->Tag()), leaf_count,
                      parent_cut_direction);
      auto child0_bitset = tag_to_bitset.at(child0_node->Tag());
      if (child0_direction) child0_bitset.flip();
      auto child1_bitset = tag_to_bitset.at(child1_node->Tag());
      if (child1_direction) child1_bitset.flip();
      bitset.CopyFrom(std::min(child0_bitset, child1_bitset), 2 * leaf_count,
                      false);
      subsplit_support.increment(bitset, count);
    });
  }
  return subsplit_support;
}

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Build") {
  Driver driver;

  auto t = driver.ParseString("((0,1),(2,(3,4)));");
  auto m = TagBitsetMapOf(t);

  std::cout << t->Newick() << std::endl;
  // TODO Add an actual test.
  PrintTagBitsetMap(m);
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BUILD_HPP_
