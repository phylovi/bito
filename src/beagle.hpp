// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BEAGLE_HPP_
#define SRC_BEAGLE_HPP_

#include <string>
#include "libhmsbeagle/beagle.h"
#include "typedefs.hpp"

namespace beagle {

// DNA assumption here.
CharIntMap GetSymbolTable() {
  CharIntMap table({{'A', 0},
                    {'C', 1},
                    {'G', 2},
                    {'T', 3},
                    {'a', 0},
                    {'c', 1},
                    {'g', 2},
                    {'t', 3},
                    {'-', 4}});
  return table;
  }

  SymbolVector SymbolVectorOf(const std::string &str,
                              const CharIntMap &symbol_table) {
    SymbolVector v(str.size());
    for (size_t i = 0; i < str.size(); i++) {
      auto search = symbol_table.find(str[i]);
      if (search != symbol_table.end()) {
        v[i] = search->second;
    } else {
      std::cerr << "Symbol '" << search->first << "' not known.\n";
      abort();
    }
    }
    return v;
  }

  int CreateInstance(int tip_count, int alignment_length,
                     BeagleInstanceDetails *return_info) {
    // Number of partial buffers to create (input) -- internal node count
    int partials_buffer_count = tip_count - 1;
    // Number of compact state representation buffers to create -- for use with
    // setTipStates (input) */
    int compact_buffer_count = tip_count;
    // DNA assumption here.
    int state_count = 4;
    // Number of site patterns to be handled by the instance (input) -- not
    // compressed in this case
    int pattern_count = alignment_length;
    // Number of eigen-decomposition buffers to allocate (input)
    int eigen_buffer_count = 1;
    // Number of transition matrix buffers (input) -- one per edge
    int matrix_buffer_count = 2 * tip_count - 1;
    // Number of rate categories
    int category_count = 1;
    // Number of scaling buffers -- can be zero if scaling is not needed
    int scale_buffer_count = 0;
    // List of potential resources on which this instance is allowed (input,
    // NULL implies no restriction
    int *allowed_resources = nullptr;
    // Length of resourceList list (input) -- not needed to use the default
    // hardware config
    int resource_count = 0;
    // Bit-flags indicating preferred implementation charactertistics, see
    // BeagleFlags (input)
    long preference_flags = 0;
    // Bit-flags indicating required implementation characteristics, see
    // BeagleFlags (input)
    int requirement_flags = 0;

    return beagleCreateInstance(
        tip_count, partials_buffer_count, compact_buffer_count, state_count,
        pattern_count, eigen_buffer_count, matrix_buffer_count, category_count,
        scale_buffer_count, allowed_resources, resource_count, preference_flags,
        requirement_flags, return_info);
  }
  }  // namespace beagle

#ifdef DOCTEST_LIBRARY_INCLUDED
  TEST_CASE("Beagle") {
    CharIntMap symbol_table = beagle::GetSymbolTable();
    SymbolVector symbol_vector =
        beagle::SymbolVectorOf("-tgcaTGCA", symbol_table);
    SymbolVector correct_symbol_vector = {4, 3, 2, 1, 0, 3, 2, 1, 0};
    CHECK_EQ(symbol_vector, correct_symbol_vector);
  }
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_BEAGLE_HPP_
