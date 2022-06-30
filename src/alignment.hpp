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
  // Return map of taxon names to sequence alignments.
  StringStringMap Data() const { return data_; }
  // Number of taxon sequences in data map.
  size_t SequenceCount() const { return data_.size(); }
  // The length of the sequence alignments.
  size_t Length() const;
  // Compare if alignments have same name and sequence data.
  bool operator==(const Alignment& other) const { return data_ == other.Data(); }
  // Is the alignment non-empty and do all sequences have the same length?
  bool IsValid() const;
  // Get alignment sequence by taxon name.
  const std::string& at(const std::string& taxon) const;
  // Load fasta file into Alignment.
  static Alignment ReadFasta(const std::string& fname);
  // Create a new alignment
  Alignment ExtractSingleColumnAlignment(size_t which_column) const;

  static Alignment HelloAlignment() {
    return Alignment({{"mars", "CCGAG-AGCAGCAATGGAT-GAGGCATGGCG"},
                      {"saturn", "GCGCGCAGCTGCTGTAGATGGAGGCATGACG"},
                      {"jupiter", "GCGCGCAGCAGCTGTGGATGGAAGGATGACG"}});
  }

 private:
  // - Map of alignments: [ taxon name -> alignment sequence ]
  StringStringMap data_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Alignment") {
  auto alignment = Alignment::ReadFasta("data/hello.fasta");
  CHECK_EQ(alignment, Alignment::HelloAlignment());
  CHECK(alignment.IsValid());
  CHECK_THROWS(alignment.ExtractSingleColumnAlignment(31));
  auto first_col_expected =
      Alignment({{"mars", "C"}, {"saturn", "G"}, {"jupiter", "G"}});
  CHECK_EQ(alignment.ExtractSingleColumnAlignment(0), first_col_expected);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
