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

using BitsetVector = std::vector<Bitset>;
using SizeBitsetMap = std::unordered_map<size_t, Bitset>;
using BitsetSizeMap = std::unordered_map<Bitset, size_t>;
using BitsetSizePairMap = std::unordered_map<Bitset, std::pair<size_t, size_t>>;
using BitsetSizeDict = DefaultDict<Bitset, size_t>;
using IndexerRepresentation = std::pair<SizeVector, SizeVectorVector>;
using PCSSDict = std::unordered_map<Bitset, DefaultDict<Bitset, size_t>>;
using PCSSIndexVector = std::vector<size_t>;

using StringSizePairMap = std::unordered_map<std::string, std::pair<size_t, size_t>>;
using SizeStringMap = std::unordered_map<size_t, std::string>;
using StringPCSSMap =
    std::unordered_map<std::string, std::unordered_map<std::string, size_t>>;

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

namespace SBNMaps {
SizeBitsetMap IdIdSetMapOf(Node::NodePtr topology);

BitsetSizeDict RootsplitCounterOf(const Node::TopologyCounter& topologies);
PCSSDict PCSSCounterOf(const Node::TopologyCounter& topologies);
// This function returns a vector indexed by the edges of the tree and
// containing the corresponding split index as indexed by the indexer.
SizeVector SplitIndicesOf(const BitsetSizeMap& indexer, const Node::NodePtr& topology);
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

StringPCSSMap StringPCSSMapOf(PCSSDict d);

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
