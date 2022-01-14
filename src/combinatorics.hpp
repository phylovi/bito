// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <cmath>
#include <cstddef>

namespace Combinatorics {

// The number of topologies on the given number of taxa.
double TopologyCount(size_t taxon_count);
// The log of the number of topologies on the given number of taxa.
double LogTreeCount(size_t taxon_count);

// Define the child subsplit count ratio for (n0, n1) as the number of topologies with
// n0 taxa times the number of topologies with n1 taxa divided by the number of
// topologies with n0+n1 taxa. This is a prior probability for a subsplit with (n0, n1)
// taxa conditioned on it resolving a subsplit on n0+n1 taxa under the uniform
// distribution on topologies.
//
// Naive version:
double LogChildSubsplitCountRatioNaive(size_t child0_taxon_count,
                                       size_t child1_taxon_count);
// Non-naive version:
double LogChildSubsplitCountRatio(size_t child0_taxon_count, size_t child1_taxon_count);

}  // namespace Combinatorics

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
