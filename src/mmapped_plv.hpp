// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// RAII class for partial likelihood vectors that are mmapped to disk.
//
// This class is to allocate a very large partial likelihood vector in virtual memory
// and then cut it up (via Subdivide) into a vector of partial likelihood vectors.

#pragma once

#include "eigen_sugar.hpp"
#include "mmapped_matrix.hpp"

using NucleotidePLV = Eigen::Matrix<double, 4, Eigen::Dynamic, Eigen::ColMajor>;
using NucleotidePLVRef = Eigen::Ref<NucleotidePLV>;
using NucleotidePLVRefVector = std::vector<NucleotidePLVRef>;

class MmappedNucleotidePLV {
 public:
  constexpr static Eigen::Index base_count_ = 4;

  MmappedNucleotidePLV(const std::string &file_path, Eigen::Index total_plv_length)
      : mmapped_matrix_(file_path, base_count_, total_plv_length){};

  void Resize(Eigen::Index total_plv_length) {
    mmapped_matrix_.ResizeMMap(base_count_, total_plv_length);
  }

  NucleotidePLVRefVector Subdivide(size_t into_count) {
    Assert(into_count > 0, "into_count is zero in MmappedNucleotidePLV::Subdivide.");
    auto entire_plv = mmapped_matrix_.Get();
    const auto total_plv_length = entire_plv.cols();
    Assert(total_plv_length % into_count == 0,
           "into_count isn't a multiple of total PLV length in "
           "MmappedNucleotidePLV::Subdivide.");
    const size_t block_length = total_plv_length / into_count;
    NucleotidePLVRefVector sub_plvs;
    sub_plvs.reserve(into_count);
    for (size_t idx = 0; idx < into_count; ++idx) {
      sub_plvs.push_back(
          entire_plv.block(0, idx * block_length, base_count_, block_length));
    }
    return sub_plvs;
  }

  size_t ByteCount() const { return mmapped_matrix_.ByteCount(); }

 private:
  MmappedMatrix<NucleotidePLV> mmapped_matrix_;
};
