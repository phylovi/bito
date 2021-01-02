// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_COUNT_HPP_
#define SRC_TREE_COUNT_HPP_

#include <cmath>
#include <cstddef>

namespace TreeCount {
double TreeCount(size_t taxon_count);
double LogTreeCount(size_t taxon_count);
// TODO(e) docs
double LogChildSubsplitCountRatioNaive(size_t child0_taxon_count,
                                       size_t child1_taxon_count);
double LogChildSubsplitCountRatio(size_t child0_taxon_count, size_t child1_taxon_count);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TreeCount") {
  CHECK_EQ(TreeCount::TreeCount(1), 1.);
  CHECK_EQ(TreeCount::TreeCount(2), 1.);
  CHECK_EQ(TreeCount::TreeCount(3), 3.);
  CHECK_EQ(TreeCount::TreeCount(4), 15.);
  CHECK_EQ(TreeCount::TreeCount(5), 105.);
  CHECK_EQ(TreeCount::TreeCount(6), 945.);
  CHECK_EQ(TreeCount::TreeCount(7), 10395.);

  for (size_t taxon_count = 1; taxon_count < 20; taxon_count++) {
    CHECK_LT(fabs(TreeCount::LogTreeCount(taxon_count) -
                  std::log(TreeCount::TreeCount(taxon_count))),
             1e-10);
  }

  // TODO(e) add a non-naive version and a test for it.
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_COUNT_HPP_
