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
typedef std::unordered_map<Bitset, int> BitsetIndexer;
typedef DefaultDict<Bitset, uint32_t> BitsetUInt32Dict;
typedef std::unordered_map<uint32_t, Bitset> UInt32BitsetMap;
typedef std::unordered_map<Bitset, uint32_t> BitsetUInt32Map;
typedef std::unordered_map<Bitset, std::pair<uint32_t, uint32_t>>
    BitsetUInt32PairMap;

typedef std::unordered_map<Bitset, DefaultDict<Bitset, uint32_t>> PCSSDict;

TagBitsetMap TagLeafSetMapOf(Node::NodePtr topology);
TagBitsetMap TagIndexSetMapOf(Node::NodePtr topology);
void PrintTagBitsetMap(TagBitsetMap map);

BitsetUInt32Dict RootsplitCounterOf(const Node::TopologyCounter& topologies);
PCSSDict PCSSCounterOf(const Node::TopologyCounter& topologies);

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Build") {
  auto topology3 = Node::ExampleTopologies()[3];

  for (auto topology : Node::ExampleTopologies()) {
    std::cout << topology->Newick([](const Node* node) {
      return std::to_string(node->Index());
    }) << std::endl;
  }

  std::cout << topology3->Newick([](const Node* node) {
    return std::to_string(node->Index());
  }) << std::endl;

  PrintTagBitsetMap(TagIndexSetMapOf(topology3));

  // 3_3 0111110
  // 3_4 1111111
  // 1_1 0100000
  // 0_1 1000000
  // 2_1 0010000
  // 3_1 0001000
  // 3_2 0011100

  // Tests comparing to vbpi appear in Python test code.
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_BUILD_HPP_
