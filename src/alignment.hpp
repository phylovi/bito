// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ALIGNMENT_HPP_
#define SRC_ALIGNMENT_HPP_

#include <string>
#include "sugar.hpp"

class Alignment {
 public:
  Alignment() = default;
  explicit Alignment(StringStringMap data) : data_(std::move(data)) {}

  StringStringMap Data() const { return data_; }
  size_t SequenceCount() const { return data_.size(); }
  size_t Length() const;
  bool operator==(const Alignment& other) const {
    return data_ == other.Data();
  }

  // Is the alignment non-empty and do all sequences have the same length?
  bool IsValid() const;
  std::string at(const std::string& taxon) const;

  static Alignment ReadFasta(const std::string& fname);

 private:
  StringStringMap data_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Alignment") {
  auto alignment = Alignment::ReadFasta("data/hello.fasta");
  Alignment correct({{"mars", "CCGAG-AGCAGCAATGGAT-GAGGCATGGCG"},
                     {"saturn", "GCGCGCAGCTGCTGTAGATGGAGGCATGACG"},
                     {"jupiter", "GCGCGCAGCAGCTGTGGATGGAAGGATGACG"}});
  CHECK_EQ(correct, alignment);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ALIGNMENT_HPP_
