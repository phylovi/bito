// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "mersenne_twister.hpp"

std::random_device MersenneTwister::random_device_;
std::mt19937 MersenneTwister::random_generator_(MersenneTwister::random_device_());
