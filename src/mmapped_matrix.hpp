// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// RAII class for Eigen matrices that are mmapped to disk.
//
// https://en.wikipedia.org/wiki/Memory-mapped_file
// Motivating post: https://stackoverflow.com/a/43301909/467327
// Intro: https://www.youtube.com/watch?v=m7E9piHcfr4
// Simple example: https://jameshfisher.com/2017/01/28/mmap-file-write/
// Most complete example: https://gist.github.com/marcetcheverry/991042

#ifndef SRC_MMAPPED_MATRIX_HPP_
#define SRC_MMAPPED_MATRIX_HPP_

#include <errno.h>
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
  MmappedMatrix(const std::string &file_path, Eigen::Index rows, Eigen::Index cols)
      : rows_(rows), cols_(cols), byte_count_(rows * cols * sizeof(Scalar)) {
    file_descriptor_ = open(
        file_path.c_str(),
        O_RDWR | O_CREAT,  // Open for reading and writing; create if it doesn't exit.
        S_IRUSR | S_IWUSR  // Make the file readable and writable by the user.
    );
    if (file_descriptor_ == -1) {
      Failwith("MmappedMatrix could not create a file at " + file_path);
    }
    // Resizes file so it's just right for our vector.
    auto ftruncate_status = ftruncate(file_descriptor_, byte_count_);
    if (ftruncate_status != 0) {
      Failwith("MmappedMatrix could not resize the file at " + file_path);
    }
    mmapped_memory_ = static_cast<Scalar *>(mmap(  //
        NULL,                    // This address is ignored as we are using MAP_SHARED.
        byte_count_,             // Size of map.
        PROT_READ | PROT_WRITE,  // We want to read and write.
        MAP_SHARED,              // We need MAP_SHARED to actually write to memory.
        file_descriptor_,        // File descriptor.
        0                        // Offset.
        ));
    if (mmapped_memory_ == MAP_FAILED) {
      throw std::system_error(errno, std::system_category(), "mmap");
    }
  }

  ~MmappedMatrix() {
    auto CheckStatus = [](int status, std::string name) {
      if (status != 0) {
        std::cout << "Warning: " << name
                  << " did not succeed in MmappedMatrix: " << strerror(errno)
                  << std::endl;
      }
    };
    // Synchronize memory with physical storage.
    auto msync_status = msync(mmapped_memory_, byte_count_, MS_SYNC);
    CheckStatus(msync_status, "msync");
    // Unmap memory mapped with mmap.
    auto munmap_status = munmap(mmapped_memory_, byte_count_);
    CheckStatus(munmap_status, "munmap");
    auto close_status = close(file_descriptor_);
    CheckStatus(close_status, "close");
  }

  MmappedMatrix(const MmappedMatrix &) = delete;
  MmappedMatrix(const MmappedMatrix &&) = delete;
  MmappedMatrix &operator=(const MmappedMatrix &) = delete;
  MmappedMatrix &operator=(const MmappedMatrix &&) = delete;

  Eigen::Map<EigenDenseMatrixBaseT> Get() {
    return Eigen::Map<EigenDenseMatrixBaseT>(mmapped_memory_, rows_, cols_);
  }

  size_t ByteCount() const { return byte_count_; }

 private:
  Eigen::Index rows_;
  Eigen::Index cols_;
  size_t byte_count_;
  int file_descriptor_;
  Scalar *mmapped_memory_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("MmappedMatrix") {
  Eigen::Index rows = 4;
  Eigen::Index cols = 5;
  using MmappedMatrixXd = MmappedMatrix<Eigen::MatrixXd>;
  {
    MmappedMatrixXd mmapped_matrix("_ignore/mmapped_matrix.data", rows, cols);
    mmapped_matrix.Get()(rows - 1, cols - 1) = 5.;
  }  // End of scope, so our mmap is destroyed and file written.
  MmappedMatrixXd mmapped_matrix("_ignore/mmapped_matrix.data", rows, cols);
  CHECK_EQ(mmapped_matrix.Get()(rows - 1, cols - 1), 5.);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_MMAPPED_MATRIX_HPP_
