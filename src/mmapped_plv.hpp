// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// RAII class for partial likelihood vectors that are mmapped to disk.

#ifndef SRC_MMAPPED_PLV_HPP_
#define SRC_MMAPPED_PLV_HPP_

#include "eigen_sugar.hpp"
#include "mmapped_matrix.hpp"

using NucleotidePLV = Eigen::Matrix<double, 4, Eigen::Dynamic, Eigen::ColMajor>;
using NucleotidePLVMap = Eigen::Map<NucleotidePLV>;
using NucleotidePLVRef = Eigen::Ref<NucleotidePLV>;
using NucleotidePLVRefVector = std::vector<NucleotidePLVRef>;

class MmappedNucleotidePLV {
 public:
  MmappedNucleotidePLV(){};
  MmappedNucleotidePLV(std::string file_path, Eigen::Index total_plv_length)
      : mmapped_matrix_(file_path, 4, total_plv_length),
        total_plv_length_(total_plv_length){};

  MmappedNucleotidePLV(const MmappedNucleotidePLV &) = delete;
  MmappedNucleotidePLV(const MmappedNucleotidePLV &&) = delete;
  MmappedNucleotidePLV &operator=(const MmappedNucleotidePLV &) = delete;
  MmappedNucleotidePLV &operator=(const MmappedNucleotidePLV &&) = delete;

  NucleotidePLVRefVector Subdivide(size_t into_count) {
    Assert(total_plv_length_ % into_count == 0,
           "into_count isn't a multiple of total PLV length in "
           "MmappedNucleotidePLV::Subdivide.");
    size_t block_length = total_plv_length_ / into_count;
    NucleotidePLVRefVector sub_plvs;
    sub_plvs.reserve(into_count);
    auto entire_plv = mmapped_matrix_.Get();
    for (size_t idx = 0; idx < into_count; ++idx) {
      sub_plvs.push_back(entire_plv.block(0, idx * block_length, 4, block_length));
    }
    return sub_plvs;
  }

 private:
  MmappedMatrix<NucleotidePLV> mmapped_matrix_;
  Eigen::Index total_plv_length_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("MmappedNucleotidePLV") {
  MmappedNucleotidePLV mmapped_plv("_ignore/mmapped_plv.data", 10);
  auto plvs = mmapped_plv.Subdivide(2);
  for (const auto &plv : plvs) {
    CHECK(plv.rows() == 4);
    CHECK(plv.cols() == 5);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_MMAPPED_PLV_HPP_
