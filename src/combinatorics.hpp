// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_COMBINATORICS_HPP_
#define SRC_COMBINATORICS_HPP_

#include <cmath>
#include <cstddef>
#include <vector>

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

// Take the cartesian product of a vector of vectors.
// Simple. Not especially efficient.
// https://stackoverflow.com/a/17050528/467327
template <typename T>
std::vector<std::vector<T>> CartesianProduct(
    const std::vector<std::vector<T>>& input_vector_vector) {
  if (input_vector_vector.empty()) {
    return {};
  }  // else
  std::vector<std::vector<T>> result = {{}};
  for (const auto& current_vector : input_vector_vector) {
    std::vector<std::vector<T>> updated_result;
    for (const auto& partial_element_of_result : result) {
      for (const auto item : current_vector) {
        updated_result.push_back(partial_element_of_result);
        updated_result.back().push_back(item);
      }
    }
    result = std::move(updated_result);
  }
  return result;
}

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

  std::vector<std::vector<int>> v = {{1, 2}, {3, 4, 5}};
  std::vector<std::vector<int>> v_result = {{1, 3}, {1, 4}, {1, 5},
                                            {2, 3}, {2, 4}, {2, 5}};
  CHECK_EQ(v_result, Combinatorics::CartesianProduct(v));
  std::vector<std::vector<int>> empty_case = {};
  CHECK_EQ(empty_case, Combinatorics::CartesianProduct(empty_case));
  std::vector<std::vector<int>> w = {{1, 2}};
  std::vector<std::vector<int>> w_result = {{1}, {2}};
  CHECK_EQ(w_result, Combinatorics::CartesianProduct(w));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_COMBINATORICS_HPP_
