// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Interface to the zlib compression library

#pragma once

#include <zlib.h>

#include <atomic>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <streambuf>

namespace zlib {

// The result of an inflate operation. Contains the counts of consumed and
// produced bytes.
struct Result {
  enum class [[nodiscard]] Code : decltype(Z_OK){
      ok = Z_OK,
      stream_end = Z_STREAM_END,
      need_dict = Z_NEED_DICT,
  };

  Code code;
  size_t in_count;
  size_t out_count;
};

// Thrown when zlib generates an error internally
struct Exception : public std::runtime_error {
  enum class Code : decltype(Z_OK) {
    sys = Z_ERRNO,
    stream = Z_STREAM_ERROR,
    data = Z_DATA_ERROR,
    mem = Z_MEM_ERROR,
    buf = Z_BUF_ERROR,
    version = Z_VERSION_ERROR,
  };

  Exception(Code c, const char* message) : std::runtime_error{message}, code{c} {}

  explicit Exception(Code c) : std::runtime_error{""}, code{c} {}

  const Code code;
};

// Flush mode; see zlib documentation. We're only using partial for maximum
// compatability.
enum class Flush : decltype(Z_OK) {
  no = Z_NO_FLUSH,
  partial = Z_PARTIAL_FLUSH,
  sync = Z_SYNC_FLUSH,
  full = Z_FULL_FLUSH,
  finish = Z_FINISH,
  block = Z_BLOCK,
  trees = Z_TREES,
};

namespace detail {

// Wrapper for zlib calls. Will throw if unsuccessful
constexpr Result::Code call_zlib(int zlib_code, const z_stream& impl);

}  // namespace detail

// Wrapper around zlib inflate
class ZStream {
 public:
  ZStream();
  ~ZStream() noexcept;

  // Manually release the zlib resources, and throw if an error occurs. This
  // will be performed automatically by the destructor, with exceptions
  void Close();

  // Performs a round of decompression. The result object contains
  // status code, count of consumed bytes and count of produced bytes
  Result Inflate(Flush mode, const unsigned char* in, size_t in_size,
                 unsigned char* out, size_t out_size);

 private:
  ::z_stream impl_ = {};
  std::atomic_flag closed_ = ATOMIC_FLAG_INIT;
};

// IO streams interface for zlib. Performs buffering on both compressed and
// decompressed sides. Takes input from a std::istream.
class ZStringBuf : public std::stringbuf {
  using base = std::stringbuf;

 public:
  // Takes an input stream for compressed side, and buffer sizes for the
  // compressed and decompressed sides.
  ZStringBuf(const std::istream& in, size_t in_buf_size, size_t out_buf_size);

  virtual ~ZStringBuf() override;

 protected:
  virtual int_type underflow() override;
  virtual int_type uflow() override;
  virtual std::streamsize xsgetn(char_type* s, std::streamsize count) override;

 private:
  void ensure_avail(std::streamsize count);

  std::streambuf& in_;
  std::unique_ptr<char[]> in_buf_;
  const std::streamsize in_buf_size_;
  std::unique_ptr<char[]> out_buf_;
  const std::streamsize out_buf_size_;
  ZStream inflate_;
};

}  // namespace zlib
