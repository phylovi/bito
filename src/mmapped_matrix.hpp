// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// RAII class for Eigen matrices that are mmapped to disk.
//
// https://en.wikipedia.org/wiki/Memory-mapped_file
// Motivating post: https://stackoverflow.com/a/43301909/467327
// Intro: https://www.youtube.com/watch?v=m7E9piHcfr4
// Simple example: https://jameshfisher.com/2017/01/28/mmap-file-write/
// Most complete example: https://gist.github.com/marcetcheverry/991042

#pragma once

#include <errno.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <iostream>

#include "eigen_sugar.hpp"
#include "sugar.hpp"

template <typename EigenDenseMatrixBaseT>
class MmappedMatrix {
  using Scalar = typename Eigen::DenseBase<EigenDenseMatrixBaseT>::Scalar;

 public:
  MmappedMatrix(const std::string &file_path, Eigen::Index rows, Eigen::Index cols)
      : rows_(rows),
        cols_(cols),
        byte_count_(rows * cols * sizeof(Scalar)),
        file_path_(file_path) {
    file_descriptor_ = open(
        file_path.c_str(),
        O_RDWR | O_CREAT,  // Open for reading and writing; create if it doesn't exit.
        S_IRUSR | S_IWUSR  // Make the file readable and writable by the user.
    );
    if (file_descriptor_ == -1) {
      Failwith("MmappedMatrix could not create a file at " + file_path_);
    }
    InitializeMMap(rows, cols);
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

  // Intialize virtual memory map.
  void InitializeMMap(Eigen::Index rows, Eigen::Index cols) {
    // Resizes file so it's just right for our vector.
    rows_ = rows;
    cols_ = cols;
    byte_count_ = rows * cols * sizeof(Scalar);
    auto ftruncate_status = ftruncate(file_descriptor_, byte_count_);

    if (ftruncate_status != 0) {
      Failwith("MmappedMatrix could not intialize the file at " + file_path_);
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

  // Resize virtual memory map. Optional reporting.
  void ResizeMMap(Eigen::Index rows, Eigen::Index cols, const bool quiet = true) {
    std::stringstream dev_null;
    auto &our_ostream = quiet ? dev_null : std::cout;
    // Resizes file so it's just right for our vector.
    rows_ = rows;
    cols_ = cols;
    const size_t old_byte_count = byte_count_;
    byte_count_ = rows * cols * sizeof(Scalar);
    if (byte_count_ == old_byte_count) {
      return;
    }
    auto ftruncate_status = ftruncate(file_descriptor_, byte_count_);
    void *old_mmapped_memory = static_cast<void *>(mmapped_memory_);

    if (ftruncate_status != 0) {
      Failwith("MmappedMatrix could not resize the file at " + file_path_);
    }

// OSX mman and mman-win32 do not implement mremap or MREMAP_MAYMOVE.
#ifndef MREMAP_MAYMOVE
    if (munmap(mmapped_memory_, old_byte_count) == -1) {
      throw std::system_error(errno, std::system_category(), "mremap");
    }

    mmapped_memory_ = static_cast<Scalar *>(mmap(  //
        NULL,                    // This address is ignored as we are using MAP_SHARED.
        byte_count_,             // Size of map.
        PROT_READ | PROT_WRITE,  // We want to read and write.
        MAP_SHARED,              // We need MAP_SHARED to actually write to memory.
        file_descriptor_,        // File descriptor.
        0                        // Offset.
        ));
#else
    mmapped_memory_ = static_cast<Scalar *>(mremap(  //
        static_cast<void *>(mmapped_memory_),        // old address
        old_byte_count,                              // old size
        byte_count_,                                 // new size
        MREMAP_MAYMOVE  // will move virtual memory if necessary.
        ));
#endif
    if (mmapped_memory_ == MAP_FAILED) {
      throw std::system_error(errno, std::system_category(), "mremap");
    }

    // Report if remapping occurred.
    bool is_remapping_moved =
        (static_cast<void *>(mmapped_memory_) != old_mmapped_memory);
    our_ostream << "REMAPPING OCCURRED: " << old_byte_count << " bytes -> "
                << byte_count_ << " bytes " << std::endl;
    our_ostream << "REMAPPING DID " << (is_remapping_moved ? "" : "*NOT* ")
                << "MOVE MEMORY ADDRESS: " << old_mmapped_memory << " "
                << mmapped_memory_ << std::endl;
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
  std::string file_path_;
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
