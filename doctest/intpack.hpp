#pragma once

#include "../src/intpack.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
inline void TestPacking(uint32_t a, uint32_t b) {
  auto p = PackInts(a, b);
  CHECK_EQ(UnpackFirstInt(p), a);
  CHECK_EQ(UnpackSecondInt(p), b);
}

TEST_CASE("intpack") {
  TestPacking(3, 4);
  TestPacking(UINT32_MAX, 4);
  TestPacking(3, UINT32_MAX);
  TestPacking(UINT32_MAX - 1, UINT32_MAX);

  // The ints are packed such that the first int takes priority in sorting.
  CHECK_LT(PackInts(0, 4), PackInts(1, 0));
}
#endif  // DOCTEST_LIBRARY_INCLUDED
