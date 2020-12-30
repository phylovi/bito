// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TREE_COUNT_HPP_
#define SRC_TREE_COUNT_HPP_

#include <cmath>
#include <cstddef>

namespace TreeCount {
double TreeCount(size_t taxon_count);
double LogTreeCount(size_t taxon_count);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TreeCount") {
  for (size_t taxon_count = 1; taxon_count < 20; taxon_count++) {
    CHECK_LT(fabs(TreeCount::LogTreeCount(taxon_count) -
                  std::log(TreeCount::TreeCount(taxon_count))),
             1e-10);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TREE_COUNT_HPP_
