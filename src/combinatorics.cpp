// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "combinatorics.hpp"

double Combinatorics::TopologyCount(size_t taxon_count) {
  double result = 1.;

  for (size_t i = 2; i <= taxon_count; i++) {
    result *= 2. * static_cast<double>(i) - 3.;
  }
  return result;
}

double Combinatorics::LogTreeCount(size_t taxon_count) {
  double result = 0.;

  for (size_t i = 2; i <= taxon_count; i++) {
    result += std::log(2. * static_cast<double>(i) - 3.);
  }
  return result;
}

double Combinatorics::LogChildSubsplitCountRatioNaive(size_t child0_taxon_count,
                                                      size_t child1_taxon_count) {
  return LogTreeCount(child0_taxon_count) + LogTreeCount(child1_taxon_count) -
         LogTreeCount(child0_taxon_count + child1_taxon_count);
}

double Combinatorics::LogChildSubsplitCountRatio(size_t child0_taxon_count,
                                                 size_t child1_taxon_count) {
  size_t total_count = child0_taxon_count + child1_taxon_count;
  double total_without_child1 = 0.;

  for (size_t i = child1_taxon_count + 1; i <= total_count; i++) {
    total_without_child1 += std::log(2. * static_cast<double>(i) - 3.);
  }
  return LogTreeCount(child0_taxon_count) - total_without_child1;
}
