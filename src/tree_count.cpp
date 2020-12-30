// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "tree_count.hpp"

double TreeCount::TreeCount(size_t taxon_count) {
  double result = 1.;

  for (size_t i = 2; i <= taxon_count; i++) {
    result *= 2. * static_cast<double>(i) - 3.;
  }
  return result;
}

double TreeCount::LogTreeCount(size_t taxon_count) {
  double result = 0.;

  for (size_t i = 2; i <= taxon_count; i++) {
    result += std::log(2. * static_cast<double>(i) - 3.);
  }
  return result;
}
