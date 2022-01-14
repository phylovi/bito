// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "mersenne_twister.hpp"

std::random_device MersenneTwister::random_device_;
std::mt19937 MersenneTwister::random_generator_(MersenneTwister::random_device_());
