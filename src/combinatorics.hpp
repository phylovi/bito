// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_COMBINATORICS_HPP_
#define SRC_COMBINATORICS_HPP_

#include <cmath>
#include <cstddef>

namespace Combinatorics {
double TopologyCount(size_t taxon_count);
double LogTreeCount(size_t taxon_count);
// TODO(e) docs
double LogChildSubsplitCountRatioNaive(size_t child0_taxon_count,
                                       size_t child1_taxon_count);
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

  // TODO(e) add a non-naive version and a test for it.
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_COMBINATORICS_HPP_
