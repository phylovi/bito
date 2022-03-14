// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "zlib_stream.hpp"

namespace zlib {

namespace detail {

constexpr Result::Code call_zlib(int zlib_code, const z_stream& impl) {
  if (zlib_code < 0) {
    const auto code = static_cast<typename Exception::Code>(zlib_code);
    if (impl.msg) {
      throw Exception(code, impl.msg);
    } else {
      throw Exception(code);
    }
  }
  return static_cast<Result::Code>(zlib_code);
}

}  // namespace detail

ZStream::ZStream() {
  impl_.next_in = nullptr;
  impl_.avail_in = 0;
  impl_.zalloc = nullptr;
  impl_.zfree = nullptr;
  impl_.opaque = this;
  auto ret = detail::call_zlib(::inflateInit2(&impl_, 16 + MAX_WBITS), impl_);
  if (ret != Result::Code::ok) {
    throw std::logic_error("Unexpected result from zlib");
  }
}

ZStream::~ZStream() noexcept try { Close(); } catch (...) {
}

void ZStream::Close() {
  if (!closed_.test_and_set()) {
    auto ret = detail::call_zlib(::inflateEnd(&impl_), impl_);
    if (ret != Result::Code::ok) {
      throw std::logic_error("Unexpected result from zlib");
    }
  }
}

Result ZStream::Inflate(Flush mode, const unsigned char* in, size_t in_size,
                        unsigned char* out, size_t out_size) {
  impl_.next_in = const_cast<unsigned char*>(in);
  impl_.avail_in = in_size;
  impl_.next_out = out;
  impl_.avail_out = out_size;
  const auto ret = detail::call_zlib(
      ::inflate(&impl_, static_cast<std::underlying_type_t<Flush>>(mode)), impl_);
  return {ret, in_size - impl_.avail_in, out_size - impl_.avail_out};
}

ZStringBuf::ZStringBuf(const std::istream& in, size_t in_buf_size, size_t out_buf_size)
    : base{},
      in_{*in.rdbuf()},
      in_buf_{std::make_unique<char[]>(in_buf_size)},
      in_buf_size_{static_cast<std::streamsize>(in_buf_size)},
      out_buf_{std::make_unique<char[]>(out_buf_size)},
      out_buf_size_{static_cast<std::streamsize>(out_buf_size)},
      inflate_{} {}

ZStringBuf::~ZStringBuf() {}

ZStringBuf::int_type ZStringBuf::underflow() {
  ensure_avail(1);
  return base::underflow();
}

ZStringBuf::int_type ZStringBuf::uflow() {
  ensure_avail(1);
  return base::uflow();
}

std::streamsize ZStringBuf::xsgetn(char_type* s, std::streamsize count) {
  ensure_avail(count);
  return base::xsgetn(s, count);
}

void ZStringBuf::ensure_avail(std::streamsize count) {
  while (count > in_avail()) {
    // Reads from the compressed stream
    const size_t in_count =
        in_.sgetn(in_buf_.get(), std::min(in_buf_size_, count - in_avail()));

    size_t consumed = 0;
    while (consumed < in_count) {
      // Performs zlib decompression
      auto ret = inflate_.Inflate(
          Flush::partial, reinterpret_cast<unsigned char*>(in_buf_.get() + consumed),
          static_cast<size_t>(in_count - consumed),
          reinterpret_cast<unsigned char*>(out_buf_.get()),
          static_cast<size_t>(out_buf_size_));
      if (ret.code != Result::Code::ok) {
        if (ret.code == Result::Code::need_dict) {
          throw std::runtime_error("Preset dictionary is needed");
        }
        if (ret.code == Result::Code::stream_end) {
          break;
        }
        throw std::logic_error("Unexpected result from zlib");
      }
      consumed += ret.in_count;

      // Stores the decompressed data
      std::streamsize produced = 0;
      while (produced < static_cast<std::streamsize>(ret.out_count)) {
        const auto prod = sputn(out_buf_.get() + produced, ret.out_count - produced);
        if (prod < 1) break;
        produced += prod;
      }
    }
    if (consumed < 1) break;
  }
}

}  // namespace zlib
