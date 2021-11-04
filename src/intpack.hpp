// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <cstdint>
#include <string>

inline uint64_t PackInts(uint32_t a, uint32_t b) {
  return (uint64_t)((((uint64_t)a) << 32) + (uint64_t)b);
}

inline uint32_t UnpackFirstInt(uint64_t x) { return (uint32_t)(((uint64_t)x) >> 32); }

inline uint32_t UnpackSecondInt(uint64_t x) {
  return (uint32_t)(((uint64_t)x) & 0xffffffff);
}

inline std::string StringOfPackedInt(uint64_t x) {
  return (std::to_string(UnpackFirstInt(x)) + "_" + std::to_string(UnpackSecondInt(x)));
}

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
