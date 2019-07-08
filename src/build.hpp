// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BUILD_HPP_
#define SRC_BUILD_HPP_

#include <algorithm>
#include <cassert>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include "bitset.hpp"
#include "default_dict.hpp"
#include "driver.hpp"
#include "tree.hpp"

typedef std::unordered_map<uint64_t, Bitset> TagBitsetMap;
typedef std::unordered_set<Bitset> BitsetSet;
typedef std::unordered_map<Bitset, int> BitsetIndexer;
typedef DefaultDict<Bitset, uint32_t> BitsetUInt32Map;

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

BitsetUInt32Map RootsplitSupportOf(TreeCollection::TreePtrCounterPtr trees) {
  BitsetUInt32Map rootsplit_counter(0);
  for (auto iter = trees->begin(); iter != trees->end(); ++iter) {
    auto tree = iter->first;
    auto count = iter->second;
    auto tag_to_bitset = TagBitsetMapOf(tree->Root());
    auto Aux = [&rootsplit_counter, &tag_to_bitset, &count](const Node* n) {
      auto split = tag_to_bitset.at(n->Tag()).copy();
      split.Minorize();
      rootsplit_counter.increment(std::move(split), count);
    };
    for (auto child : tree->Root()->Children()) {
      child->PreOrder(Aux);
    }
  }
  return rootsplit_counter;
}

BitsetUInt32Map SubsplitSupportOf(TreeCollection::TreePtrCounterPtr trees) {
  BitsetUInt32Map subsplit_support(0);
  for (auto iter = trees->begin(); iter != trees->end(); ++iter) {
    auto tree = iter->first;
    auto count = iter->second;
    auto tag_to_bitset = TagBitsetMapOf(tree->Root());
    auto leaf_count = tree->LeafCount();
    // TODO(ematsen) make a more informative error message when people don't put
    // in a bifurcating tree with a trifurcation at the root.
    tree->Root()->PCSSPreOrder(
        [&subsplit_support, &tag_to_bitset, &count, &leaf_count](
            const Node* uncut_parent_node, bool uncut_parent_direction,
            const Node* cut_parent_node, bool cut_parent_direction,  //
            const Node* child0_node, bool child0_direction,          //
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

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Build") {
  Driver driver;

  auto trees = driver.ParseFile("data/many_rootings.tre")->Trees();
  auto support = SubsplitSupportOf(trees);
  // Get the support of the first tree in trees.
  auto single_tree = std::make_shared<TreeCollection::TreePtrCounter>();
  single_tree->insert(std::make_pair(trees->begin()->first, 1));
  auto single_support = SubsplitSupportOf(single_tree);
  // many_rootings has many (unrooted) rootings of the same tree.
  // Here we check to make sure that every support across the various rootings
  // is in the SBN support for the single tree.
  for (auto iter = support.begin(); iter != support.end(); ++iter) {
    CHECK(iter->first.PCSSIsValid());
    CHECK(single_support.contains(iter->first));
  }
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BUILD_HPP_
