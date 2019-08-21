// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SBN_MAPS_HPP_
#define SRC_SBN_MAPS_HPP_

#include <unordered_map>
#include <utility>
#include <vector>
#include "bitset.hpp"
#include "default_dict.hpp"
#include "driver.hpp"
#include "tree.hpp"

typedef std::vector<Bitset> BitsetVector;
typedef std::unordered_map<size_t, Bitset> SizeBitsetMap;
typedef std::unordered_map<Bitset, int> BitsetIndexer;
typedef DefaultDict<Bitset, uint32_t> BitsetUInt32Dict;
typedef std::unordered_map<uint32_t, Bitset> UInt32BitsetMap;
typedef std::unordered_map<Bitset, uint32_t> BitsetUInt32Map;
typedef std::unordered_map<Bitset, std::pair<uint32_t, uint32_t>>
    BitsetUInt32PairMap;
typedef std::pair<SizeVector, SizeVectorVector> IndexerRepresentation;

typedef std::unordered_map<Bitset, DefaultDict<Bitset, uint32_t>> PCSSDict;

SizeBitsetMap IdIdSetMapOf(Node::NodePtr topology);

BitsetUInt32Dict RootsplitCounterOf(const Node::TopologyCounter& topologies);
PCSSDict PCSSCounterOf(const Node::TopologyCounter& topologies);
// This function gives information about the splits and PCSSs of a given
// topology with respect to the current indexing data structures.
// Specifically, it returns a pair (rootsplit_result, pcss_result).
// Each of these vectors are indexed by virtual rootings of the tree.
// rootsplit_result simply gives the indices of the rootsplits that appear for
// those various virtual rootings. pcss_result is a vector of vectors, giving
// the indices of sbn_parameters_ corresponding to PCSSs that are present in the
// given topology.
IndexerRepresentation IndexerRepresentationOf(const BitsetUInt32Map& indexer_,
                                              const Node::NodePtr& topology);

SizeVectorVector PSPRepresentationOf(const BitsetUInt32Map& indexer,
                                     const Node::NodePtr& topology);

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("SBNMaps") {
  auto topology0 = Node::ExampleTopologies()[0];

  // (0,1,(2,3)4)5;
  auto correct_id_id_set_map =
      std::unordered_map<size_t, Bitset>({{5, Bitset("111111")},
                                          {1, Bitset("010000")},
                                          {0, Bitset("100000")},
                                          {2, Bitset("001000")},
                                          {3, Bitset("000100")},
                                          {4, Bitset("001110")}});

  for (const auto& iter : IdIdSetMapOf(topology0)) {
    CHECK_EQ(correct_id_id_set_map.at(iter.first), iter.second);
  }

  // Tests comparing to vbpi appear in Python test code.
  // Tests of IndexerRepresentationOf in libsbn.hpp.
}

#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_SBN_MAPS_HPP_
