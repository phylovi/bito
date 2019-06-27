#ifndef SRC_INTPACK_HPP_
#define SRC_INTPACK_HPP_

#include <cstdint>
#include "doctest.h"

uint64_t PackInts(uint32_t a, uint32_t b) {
  return (uint64_t)((((uint64_t)a) << 32) + (uint64_t)b);
}

uint32_t UnpackFirstInt(uint64_t x) { return (uint32_t)(((uint64_t)x) >> 32); }

uint32_t UnpackSecondInt(uint64_t x) {
  return (uint32_t)(((uint64_t)x) & 0xffffffff);
}

void TestPacking(uint32_t a, uint32_t b) {
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

#endif  // SRC_INTPACK_HPP_
