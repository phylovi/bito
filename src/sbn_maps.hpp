// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A collection of functions to handle the subsplit support and to turn trees into
// indexer representations.

#ifndef SRC_SBN_MAPS_HPP_
#define SRC_SBN_MAPS_HPP_

#include <unordered_map>
#include <utility>
#include <vector>

#include "bitset.hpp"
#include "default_dict.hpp"
#include "driver.hpp"
#include "node.hpp"

using BitsetVector = std::vector<Bitset>;
using SizeBitsetMap = std::unordered_map<size_t, Bitset>;
using BitsetSizeMap = std::unordered_map<Bitset, size_t>;
using BitsetSizePairMap = std::unordered_map<Bitset, std::pair<size_t, size_t>>;
using BitsetSizeDict = DefaultDict<Bitset, size_t>;
using RootedIndexerRepresentation = SizeVector;
using RootedIndexerRepresentationCounter =
    std::vector<std::pair<RootedIndexerRepresentation, uint32_t>>;
using UnrootedIndexerRepresentation = SizeVectorVector;
using UnrootedIndexerRepresentationCounter =
    std::vector<std::pair<UnrootedIndexerRepresentation, uint32_t>>;
using PCSSDict = std::unordered_map<Bitset, DefaultDict<Bitset, size_t>>;
using PCSSIndexVector = std::vector<size_t>;
using RootedIndexerRepresentationSizeDict =
    DefaultDict<RootedIndexerRepresentation, size_t>;

using StringSizePairMap = std::unordered_map<std::string, std::pair<size_t, size_t>>;
using SizeStringMap = std::unordered_map<size_t, std::string>;
using StringPCSSMap =
    std::unordered_map<std::string, std::unordered_map<std::string, size_t>>;
using IndexerBundle =
    std::tuple<BitsetVector, BitsetSizeMap, SizeBitsetMap, BitsetSizePairMap, size_t>;

namespace SBNMaps {
// Make a map from each Tag to the bitset representing the ids below the Tag.
SizeBitsetMap IdIdSetMapOf(const Node::NodePtr& topology);
// This function returns a vector indexed by the edges of the tree and
// containing indices of the corresponding splits as indexed by the indexer.
SizeVector SplitIndicesOf(const BitsetSizeMap& indexer, const Node::NodePtr& topology);
// Make a string version of a PCSSDict.
StringPCSSMap StringPCSSMapOf(PCSSDict d);
IndexerBundle BuildIndexerBundle(const BitsetSizeDict& rootsplit_counter,
                                 const PCSSDict& pcss_counter);
}  // namespace SBNMaps

namespace UnrootedSBNMaps {
// Make a DefaultDict mapping rootsplits to the number of times they were seen.
BitsetSizeDict RootsplitCounterOf(const Node::TopologyCounter& topologies);
// Make a PCSSDict mapping PCSSs to the number of times they were seen.
PCSSDict PCSSCounterOf(const Node::TopologyCounter& topologies);
// This function gives information about the rootsplits and PCSSs of a given
// topology with respect to the current indexing data structures.
// Specifically, it returns a vector of vectors, such that the ith vector is the indices
// of sbn_parameters_ describing the tree when it is rooted above the ith node. The
// first entry of this representation is always the index of the rootsplit. The rest are
// the indices of the PCSSs that are present in the given topology.
// NOTE: Any rootsplits or PCSSs that aren't known by the indexer are assigned
// `default_index`.
UnrootedIndexerRepresentation IndexerRepresentationOf(const BitsetSizeMap& indexer,
                                                      const Node::NodePtr& topology,
                                                      const size_t default_index);
// Turn a TopologyCounter into an IndexerRepresentationCounter.
UnrootedIndexerRepresentationCounter IndexerRepresentationCounterOf(
    const BitsetSizeMap& indexer, const Node::TopologyCounter& topology_counter,
    const size_t default_index);

}  // namespace UnrootedSBNMaps

namespace RootedSBNMaps {
// Make a DefaultDict mapping rootsplits to the number of times they were seen.
BitsetSizeDict RootsplitCounterOf(const Node::TopologyCounter& topologies);
// Make a PCSSDict mapping PCSSs to the number of times they were seen.
PCSSDict PCSSCounterOf(const Node::TopologyCounter& topologies);
// A rooted indexer representation is the indexer representation of a given rooted tree.
// That is, the first entry is the rootsplit for that rooting, and after that come the
// PCSS indices.
RootedIndexerRepresentation RootedIndexerRepresentationOf(const BitsetSizeMap& indexer,
                                                          const Node::NodePtr& topology,
                                                          const size_t default_index);
// For counting standardized (i.e. PCSS index sorted) rooted indexer representations.
void IncrementRootedIndexerRepresentationSizeDict(
    RootedIndexerRepresentationSizeDict& dict,
    const RootedIndexerRepresentation rooted_indexer_representation);
// Apply the above to every rooting in the unrooted indexer representation.
void IncrementRootedIndexerRepresentationSizeDict(
    RootedIndexerRepresentationSizeDict& dict,
    const UnrootedIndexerRepresentation& indexer_representation);

}  // namespace RootedSBNMaps

// Turn a <Key, T> map into a <std::string, T> map for any Key type that has
// a ToString method.
template <class Key, class T>
std::unordered_map<std::string, T> StringifyMap(std::unordered_map<Key, T> m) {
  std::unordered_map<std::string, T> m_str;
  for (const auto& iter : m) {
    m_str[iter.first.ToString()] = iter.second;
  }
  return m_str;
}

// Hash for vectors of size_t.
// https://www.boost.org/doc/libs/1_35_0/doc/html/boost/hash_combine_id241013.html
namespace std {
template <>
struct hash<SizeVector> {
  size_t operator()(const SizeVector& values) const {
    int hash = values[0];
    for (size_t i = 1; i < values.size(); i++) {
      hash ^= values[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};
}  // namespace std

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
  // Tests of IndexerRepresentationOf in unrooted_sbn_instance.hpp.
}

#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_SBN_MAPS_HPP_
