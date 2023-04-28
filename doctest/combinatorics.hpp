#pragma once

#include "../src/combinatorics.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Combinatorics") {
  CHECK_EQ(Combinatorics::TopologyCount(1), 1.);
  CHECK_EQ(Combinatorics::TopologyCount(2), 1.);
  CHECK_EQ(Combinatorics::TopologyCount(3), 3.);
  CHECK_EQ(Combinatorics::TopologyCount(4), 15.);
  CHECK_EQ(Combinatorics::TopologyCount(5), 105.);
  CHECK_EQ(Combinatorics::TopologyCount(6), 945.);
  CHECK_EQ(Combinatorics::TopologyCount(7), 10395.);

  for (size_t taxon_count = 1; taxon_count < 20; taxon_count++) {
    CHECK_LT(fabs(Combinatorics::LogTreeCount(taxon_count) -
                  std::log(Combinatorics::TopologyCount(taxon_count))),
             1e-10);
  }

  for (size_t child0_count = 1; child0_count < 10; child0_count++) {
    for (size_t child1_count = 1; child1_count < 10; child1_count++) {
      CHECK_LT(
          fabs(Combinatorics::LogChildSubsplitCountRatio(child0_count, child1_count) -
               Combinatorics::LogChildSubsplitCountRatioNaive(child0_count,
                                                              child1_count)),
          1e-10);
    }
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
