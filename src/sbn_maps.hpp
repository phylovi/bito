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
typedef std::unordered_map<Bitset, size_t> BitsetSizeMap;
typedef std::unordered_map<Bitset, std::pair<size_t, size_t>> BitsetSizePairMap;
typedef DefaultDict<Bitset, size_t> BitsetSizeDict;
typedef std::pair<SizeVector, SizeVectorVector> IndexerRepresentation;

typedef std::unordered_map<Bitset, DefaultDict<Bitset, size_t>> PCSSDict;

namespace SBNMaps {
SizeBitsetMap IdIdSetMapOf(Node::NodePtr topology);

BitsetSizeDict RootsplitCounterOf(const Node::TopologyCounter& topologies);
PCSSDict PCSSCounterOf(const Node::TopologyCounter& topologies);
// This function returns a vector indexed by the edges of the tree and
// containing the corresponding split index as indexed by the indexer.
SizeVector SplitIndicesOf(const BitsetSizeMap& indexer,
                          const Node::NodePtr& topology);
// This function gives information about the splits and PCSSs of a given
// topology with respect to the current indexing data structures.
// Specifically, it returns a pair (rootsplit_result, pcss_result).
// Each of these vectors are indexed by virtual rootings of the tree.
// rootsplit_result simply gives the indices of the rootsplits that appear for
// those various virtual rootings. pcss_result is a vector of vectors, giving
// the indices of sbn_parameters_ corresponding to PCSSs that are present in the
// given topology.
IndexerRepresentation IndexerRepresentationOf(const BitsetSizeMap& indexer,
                                              const Node::NodePtr& topology);

// Return ragged a vector of vectors such that the ith vector is the collection
// of branch lengths in the tree collection for the ith split.
// TODO I think this is set up wrong and we want to use a PSP indexer instead of
// an indexer.
DoubleVectorVector BranchLengthsBySplit(const BitsetSizeMap& indexer,
                                        size_t split_count,
                                        const TreeCollection& tree_collection);
}  // namespace SBNMaps

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

  for (const auto& iter : SBNMaps::IdIdSetMapOf(topology0)) {
    CHECK_EQ(correct_id_id_set_map.at(iter.first), iter.second);
  }

  // Tests comparing to vbpi appear in Python test code.
  // Tests of IndexerRepresentationOf in libsbn.hpp.
}

#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_SBN_MAPS_HPP_
