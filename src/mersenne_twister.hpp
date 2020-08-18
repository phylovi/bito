// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_MERSENNE_TWISTER_HPP_
#define SRC_MERSENNE_TWISTER_HPP_

#include <random>

class MersenneTwister {
 public:
  inline void SetSeed(uint64_t seed) { random_generator_.seed(seed); }
  inline std::mt19937 &GetGenerator() const { return random_generator_; };

 private:
  static std::random_device random_device_;
  static std::mt19937 random_generator_;
};

#endif  // SRC_MERSENNE_TWISTER_HPP_
