// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ALIGNMENT_HPP_
#define SRC_ALIGNMENT_HPP_

#include "typedefs.hpp"

// The mutable nature of this class is a little ugly.
// However, its only purpose is to sit in a libsbn Instance, which is all about
// mutable state.
class Alignment {
 public:
  Alignment() {}
  explicit Alignment(StringStringMap data) : data_(data) {}

  StringStringMap Data() const { return data_; }
  size_t SequenceCount() const { return data_.size(); }
  size_t Length() const;
  // Is the alignment non-empty and do all sequences have the same length?
  bool IsValid() const;
  std::string at(const std::string &taxon) const;
  void ReadFasta(std::string fname);

 private:
  StringStringMap data_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Alignment") {
  Alignment alignment;
  alignment.ReadFasta("data/hello.fasta");
  StringStringMap correct({{"mars", "CCGAG-AGCAGCAATGGAT-GAGGCATGGCG"},
                           {"saturn", "GCGCGCAGCTGCTGTAGATGGAGGCATGACG"},
                           {"jupiter", "GCGCGCAGCAGCTGTGGATGGAAGGATGACG"}});
  CHECK_EQ(correct, alignment.Data());
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ALIGNMENT_HPP_
