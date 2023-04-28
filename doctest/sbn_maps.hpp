#pragma once

#include "../src/sbn_maps.hpp"

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
