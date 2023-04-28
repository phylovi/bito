// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <cstdint>
#include <string>

inline uint64_t PackInts(uint32_t a, uint32_t b) {
  return (static_cast<uint64_t>(a) << 32) + static_cast<uint64_t>(b);
}

inline uint32_t UnpackFirstInt(uint64_t x) { return static_cast<uint32_t>(x >> 32); }

inline uint32_t UnpackSecondInt(uint64_t x) {
  return static_cast<uint32_t>(x & 0xffffffff);
}

inline std::string StringOfPackedInt(uint64_t x) {
  return (std::to_string(UnpackFirstInt(x)) + "_" + std::to_string(UnpackSecondInt(x)));
}
