// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// RAII class for Eigen matrices that are mmapped to disk.
//
// Motivating post: https://stackoverflow.com/a/43301909/467327
// Intro: https://www.youtube.com/watch?v=m7E9piHcfr4
// Simple example: https://jameshfisher.com/2017/01/28/mmap-file-write/
// Most complete example: https://gist.github.com/marcetcheverry/991042

#ifndef SRC_MMAPPED_MATRIX_HPP_
#define SRC_MMAPPED_MATRIX_HPP_

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <Eigen/Dense>
#include <iostream>
#include "sugar.hpp"

template <typename EigenDenseMatrixBaseT>
class MmappedMatrix {
  using Scalar = typename Eigen::DenseBase<EigenDenseMatrixBaseT>::Scalar;

 public:
  MmappedMatrix(std::string file_path, Eigen::Index rows, Eigen::Index cols)
      : rows_(rows), cols_(cols) {
    file_descriptor_ = open(
        file_path.c_str(),
        O_RDWR | O_CREAT,  // Open for reading and writing; create if it doesn't exit.
        S_IRUSR | S_IWUSR  // Make the file readable and writable by the user.
    );
    if (file_descriptor_ == -1) {
      Failwith("MmappedMatrix could not create a file at " + file_path);
    }
    mmapped_len_ = rows * cols * sizeof(Scalar);
    // Resizes fd so it's just right for our vector.
    auto ftruncate_status = ftruncate(file_descriptor_, mmapped_len_);
    if (ftruncate_status != 0) {
      Failwith("MmappedMatrix could not create a file at " + file_path);
    }
    mmapped_memory_ = (Scalar *)mmap(  //
        NULL,                    // This address is ignored as we are using MAP_SHARED.
        mmapped_len_,            // Size of map.
        PROT_READ | PROT_WRITE,  // We want to read and write.
        MAP_SHARED,              // We need MAP_SHARED to actually write to memory.
        file_descriptor_,        // File descriptor.
        0                        // Offset.
    );
  }

  ~MmappedMatrix() {
    // Synchronize memory with physical storage.
    auto msync_status = msync(mmapped_memory_, mmapped_len_, MS_SYNC);
    if (msync_status != 0) {
      std::cout << "Warning: msync did not succeed in MmappedMatrix.\n";
    }
    // Unmap memory mapped with mmap.
    auto munmap_status = munmap(mmapped_memory_, mmapped_len_);
    if (munmap_status != 0) {
      std::cout << "Warning: munmap did not succeed in MmappedMatrix.\n";
    }
    auto close_status = close(file_descriptor_);
    if (close_status != 0) {
      std::cout << "Warning: close did not succeed in MmappedMatrix.\n";
    }
  }

  MmappedMatrix(const MmappedMatrix &) = delete;
  MmappedMatrix(const MmappedMatrix &&) = delete;
  MmappedMatrix &operator=(const MmappedMatrix &) = delete;
  MmappedMatrix &operator=(const MmappedMatrix &&) = delete;

  Eigen::Map<EigenDenseMatrixBaseT> get() {
    return Eigen::Map<EigenDenseMatrixBaseT>(mmapped_memory_, rows_, cols_);
  }

 private:
  int file_descriptor_;
  Eigen::Index rows_;
  Eigen::Index cols_;
  size_t mmapped_len_;
  Scalar *mmapped_memory_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("MappedMatrix") {
  using MmappedMatrixXd = MmappedMatrix<Eigen::MatrixXd>;
  MmappedMatrixXd mmapped_matrix("_ignore/mapped_matrix.data", 4, 5);
  mmapped_matrix.get()(2, 3) = 5.;
  MmappedMatrixXd mmapped_matrix_read("_ignore/mapped_matrix.data", 4, 5);
  CHECK_EQ(mmapped_matrix_read.get()(2, 3), 5.);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_MMAPPED_MATRIX_HPP_

