// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "build.hpp"
#include <algorithm>
#include <cassert>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

// Using insert and at avoids needing to make a default constructor.
// https://stackoverflow.com/questions/17172080/insert-vs-emplace-vs-operator-in-c-map

TagBitsetMap TagBitsetMapOf(Node::NodePtr t) {
  TagBitsetMap m;
  auto leaf_count = t->LeafCount();
  t->PostOrder([&m, leaf_count](const Node* n) {
    Bitset x((size_t)leaf_count);
    if (n->IsLeaf()) {
      x.set(n->MaxLeafID());
    } else {
      // Take the union of the children below.
      for (const auto& child : n->Children()) {
        x |= m.at(child->Tag());
      }
    }
    assert(m.insert({n->Tag(), std::move(x)}).second);
  });
  return m;
}

void PrintTagBitsetMap(TagBitsetMap m) {
  for (const auto& iter : m) {
    std::cout << StringOfPackedInt(iter.first) << " " << iter.second.ToString()
              << std::endl;
  }
}

BitsetUInt32Dict RootsplitCounterOf(const Node::TopologyCounter& topologies) {
  BitsetUInt32Dict rootsplit_counter(0);
  for (const auto& iter : topologies) {
    auto topology = iter.first;
    auto count = iter.second;
    auto tag_to_bitset = TagBitsetMapOf(topology);
    auto Aux = [&rootsplit_counter, &tag_to_bitset, &count](const Node* n) {
      auto split = tag_to_bitset.at(n->Tag()).copy();
      split.Minorize();
      rootsplit_counter.increment(std::move(split), count);
    };
    for (const auto& child : topology->Children()) {
      child->PreOrder(Aux);
    }
  }
  return rootsplit_counter;
}

BitsetUInt32Dict PCSSCounterOf(const Node::TopologyCounter& topologies) {
  BitsetUInt32Dict subsplit_support(0);
  for (const auto& iter : topologies) {
    auto topology = iter.first;
    auto count = iter.second;
    auto tag_to_bitset = TagBitsetMapOf(topology);
    auto leaf_count = topology->LeafCount();
    // TODO(ematsen) make a more informative error message when people don't
    // put in a bifurcating tree with a trifurcation at the root.
    topology->PCSSPreOrder(
        [&subsplit_support, &tag_to_bitset, &count, &leaf_count](
            const Node* uncut_parent_node, bool uncut_parent_direction,
            const Node* cut_parent_node,
            bool cut_parent_direction,                       //
            const Node* child0_node, bool child0_direction,  //
            const Node* child1_node, bool child1_direction) {
          Bitset bitset(3 * leaf_count, false);
          bitset.CopyFrom(tag_to_bitset.at(uncut_parent_node->Tag()), 0,
                          uncut_parent_direction);
          bitset.CopyFrom(tag_to_bitset.at(cut_parent_node->Tag()), leaf_count,
                          cut_parent_direction);
          auto child0_bitset = tag_to_bitset.at(child0_node->Tag());
          if (child0_direction) child0_bitset.flip();
          auto child1_bitset = tag_to_bitset.at(child1_node->Tag());
          if (child1_direction) child1_bitset.flip();
          bitset.CopyFrom(std::min(child0_bitset, child1_bitset),
                          2 * leaf_count, false);
          subsplit_support.increment(std::move(bitset), count);
        });
  }
  return subsplit_support;
}
