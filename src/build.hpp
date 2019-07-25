// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BUILD_HPP_
#define SRC_BUILD_HPP_

#include <unordered_map>
#include <utility>
#include <vector>
#include "bitset.hpp"
#include "default_dict.hpp"
#include "driver.hpp"
#include "tree.hpp"

typedef std::vector<Bitset> BitsetVector;
typedef std::unordered_map<uint64_t, Bitset> TagBitsetMap;
typedef std::unordered_map<size_t, Bitset> SizeBitsetMap;
typedef std::unordered_map<Bitset, int> BitsetIndexer;
typedef DefaultDict<Bitset, uint32_t> BitsetUInt32Dict;
typedef std::unordered_map<uint32_t, Bitset> UInt32BitsetMap;
typedef std::unordered_map<Bitset, uint32_t> BitsetUInt32Map;
typedef std::unordered_map<Bitset, std::pair<uint32_t, uint32_t>>
    BitsetUInt32PairMap;
typedef std::pair<SizeVector, SizeVectorVector> IndexerRepresentation;

typedef std::unordered_map<Bitset, DefaultDict<Bitset, uint32_t>> PCSSDict;

TagBitsetMap TagLeafSetMapOf(Node::NodePtr topology);
TagBitsetMap IndexIndexSetMapOf(Node::NodePtr topology);
void PrintTagBitsetMap(TagBitsetMap map);

BitsetUInt32Dict RootsplitCounterOf(const Node::TopologyCounter& topologies);
PCSSDict PCSSCounterOf(const Node::TopologyCounter& topologies);
IndexerRepresentation IndexerRepresentationOf(const BitsetUInt32Map& indexer_,
                                              const Node::NodePtr& topology);

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Build") {
  auto topology0 = Node::ExampleTopologies()[0];

  // (0,1,(2,3)4)5;
  auto correct_index_index_set_map =
      std::unordered_map<size_t, Bitset>({{5, Bitset("111111")},
                                          {1, Bitset("010000")},
                                          {0, Bitset("100000")},
                                          {2, Bitset("001000")},
                                          {3, Bitset("000100")},
                                          {4, Bitset("001110")}});

  for (const auto& iter : IndexIndexSetMapOf(topology0)) {
    CHECK_EQ(correct_index_index_set_map.at(iter.first), iter.second);
  }

  // Tests comparing to vbpi appear in Python test code.
}

#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_BUILD_HPP_
