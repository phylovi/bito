// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <string>
#include <utility>

#include "sugar.hpp"

class Alignment {
 public:
  Alignment() = default;
  explicit Alignment(StringStringMap data) : data_(std::move(data)) {}

  StringStringMap Data() const { return data_; }
  size_t SequenceCount() const { return data_.size(); }
  size_t Length() const;
  bool operator==(const Alignment& other) const { return data_ == other.Data(); }

  // Is the alignment non-empty and do all sequences have the same length?
  bool IsValid() const;
  const std::string& at(const std::string& taxon) const;

  static Alignment ReadFasta(const std::string& fname);
  Alignment ExtractSingleColumnAlignment(size_t which_column) const;

  static Alignment HelloAlignment() {
    return Alignment({{"mars", "CCGAG-AGCAGCAATGGAT-GAGGCATGGCG"},
                      {"saturn", "GCGCGCAGCTGCTGTAGATGGAGGCATGACG"},
                      {"jupiter", "GCGCGCAGCAGCTGTGGATGGAAGGATGACG"}});
  }

 private:
  StringStringMap data_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Alignment") {
  auto alignment = Alignment::ReadFasta("data/hello.fasta");
  CHECK_EQ(alignment, Alignment::HelloAlignment());
  CHECK(alignment.IsValid());
}
#endif  // DOCTEST_LIBRARY_INCLUDED
