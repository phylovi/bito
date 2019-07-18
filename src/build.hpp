// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BUILD_HPP_
#define SRC_BUILD_HPP_

#include <unordered_map>
#include <utility>
#include "bitset.hpp"
#include "default_dict.hpp"
#include "driver.hpp"
#include "tree.hpp"

typedef std::unordered_map<uint64_t, Bitset> TagBitsetMap;
typedef std::unordered_map<Bitset, int> BitsetIndexer;
typedef DefaultDict<Bitset, uint32_t> BitsetUInt32Dict;
typedef std::unordered_map<Bitset, uint32_t> BitsetUInt32Map;
typedef std::unordered_map<Bitset, std::pair<uint32_t, uint32_t>>
    BitsetUInt32PairMap;
typedef std::unordered_map<Bitset, Bitset> BitsetBitsetMap;
typedef std::unordered_map<Bitset, std::vector<Bitset>> BitsetBitsetVectorMap;

TagBitsetMap TagBitsetMapOf(Node::NodePtr t);
void PrintTagBitsetMap(TagBitsetMap m);

BitsetUInt32Dict RootsplitCounterOf(const Node::TopologyCounter& topologies);
BitsetUInt32Dict PCSSCounterOf(const Node::TopologyCounter& topologies);

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Build") {
  Driver driver;
  const auto& trees = driver.ParseNewickFile("data/many_rootings.nwk");
  auto counter = trees->TopologyCounter();
  auto support = PCSSCounterOf(counter);
  // Get the support of the first tree in trees.
  Node::TopologyCounter single_topology;
  single_topology.insert({counter.begin()->first, 1});
  auto single_support = PCSSCounterOf(single_topology);
  // many_rootings has many (unrooted) rootings of the same tree.
  // Here we check to make sure that every support across the various rootings
  // is in the SBN support for the single tree.
  for (const auto& iter : support) {
    CHECK(iter.first.PCSSIsValid());
    CHECK(single_support.contains(iter.first));
  }
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BUILD_HPP_
