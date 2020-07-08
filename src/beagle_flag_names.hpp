// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// BeagleAccessories are collections of artifacts that we can make in constant
// time given the tree, and remain const througout any operation-gathering tree
// traversal.

#ifndef SRC_BEAGLE_FLAG_NAMES_HPP_
#define SRC_BEAGLE_FLAG_NAMES_HPP_

#include <bitset>
#include <limits>
#include <string>
#include <vector>

#include "libhmsbeagle/beagle.h"

namespace BeagleFlagNames {

constexpr int long_bit_count = std::numeric_limits<long>::digits + 1;

static const std::vector<std::string> name_vector{
    "PRECISION_SINGLE",     // 00 Single precision computation
    "PRECISION_DOUBLE",     // 01 Double precision computation
    "COMPUTATION_SYNCH",    // 02 Synchronous computation (blocking)
    "COMPUTATION_ASYNCH",   // 03 Asynchronous computation (non-blocking)
    "EIGEN_REAL",           // 04 Real eigenvalue computation
    "EIGEN_COMPLEX",        // 05 Complex eigenvalue computation
    "SCALING_MANUAL",       // 06 Manual scaling
    "SCALING_AUTO",         // 07 Auto-scaling on (deprecated)
    "SCALING_ALWAYS",       // 08 Scale at every updatePartials (deprecated)
    "SCALERS_RAW",          // 09 Save raw scalers
    "SCALERS_LOG",          // 10 Save log scalers
    "VECTOR_SSE",           // 11 SSE computation
    "VECTOR_NONE",          // 12 No vector computation
    "THREADING_OPENMP",     // 13 OpenMP threading
    "THREADING_NONE",       // 14 No threading (default)
    "PROCESSOR_CPU",        // 15 Use CPU as main processor
    "PROCESSOR_GPU",        // 16 Use GPU as main processor
    "PROCESSOR_FPGA",       // 17 Use FPGA as main processor
    "PROCESSOR_CELL",       // 18 Use Cell as main processor
    "PROCESSOR_PHI",        // 19 Use Intel Phi as main processor
    "INVEVEC_STANDARD",     // 20 Inverse eigen vectors have not been transposed
    "INVEVEC_TRANSPOSED",   // 21 Inverse eigen vectors have been transposed
    "FRAMEWORK_CUDA",       // 22 Use CUDA implementation with GPU resources
    "FRAMEWORK_OPENCL",     // 23 Use OpenCL implementation with GPU resources
    "VECTOR_AVX",           // 24 AVX computation
    "SCALING_DYNAMIC",      // 25 Manual scaling with dynamic checking (deprecated)
    "PROCESSOR_OTHER",      // 26 Use other type of processor
    "FRAMEWORK_CPU",        // 27 Use CPU implementation
    "PARALLELOPS_STREAMS",  // 28 Ops may be assigned to separate device streams
    "PARALLELOPS_GRID",     // 29 Ops may be folded into single kernel launch
    "THREADING_CPP",        // 30 C++11 threading
};

std::string OfBeagleFlags(long flags) {
  std::bitset<long_bit_count> beagle_bitset{static_cast<unsigned long long>(flags)};
  std::string set_flags;
  std::string perhaps_space;
  for (size_t bit_idx = 0; bit_idx < name_vector.size(); ++bit_idx) {
    if (beagle_bitset.test(bit_idx)) {
      set_flags += perhaps_space + name_vector.at(bit_idx);
      if (perhaps_space.empty()) {
        perhaps_space = " ";
      }
    }
  }
  return set_flags;
}

}  // namespace BeagleFlagNames

#endif  // SRC_BEAGLE_FLAG_NAMES_HPP_
