#pragma once

#include "../src/site_pattern.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("SitePattern") {
  CharIntMap symbol_table = SitePattern::GetSymbolTable();
  SymbolVector symbol_vector = SitePattern::SymbolVectorOf(symbol_table, "-tgcaTGCA?");
  SymbolVector correct_symbol_vector = {4, 3, 2, 1, 0, 3, 2, 1, 0, 4};
  CHECK_EQ(symbol_vector, correct_symbol_vector);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
