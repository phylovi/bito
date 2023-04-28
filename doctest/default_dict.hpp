#pragma once

#include "../src/default_dict.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("DefaultDict") {
  auto d = DefaultDict<int, int>(0);
  CHECK_EQ(d.at(4), 0);
  d.increment(4, 5);
  CHECK_EQ(d.at(4), 5);
  d.increment(4, 2);
  CHECK_EQ(d.at(4), 7);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
