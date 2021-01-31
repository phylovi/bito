// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A class for an alignment that has been compressed into site patterns.

#ifndef SRC_SITE_PATTERN_HPP_
#define SRC_SITE_PATTERN_HPP_

#include <string>
#include <vector>

#include "alignment.hpp"
#include "sugar.hpp"

class SitePattern {
 public:
  SitePattern() = default;
  SitePattern(const Alignment& alignment, TagStringMap tag_taxon_map)
      : alignment_(alignment), tag_taxon_map_(std::move(tag_taxon_map)) {
    patterns_.resize(alignment.SequenceCount());
    Compress();
  }

  static CharIntMap GetSymbolTable();
  static SymbolVector SymbolVectorOf(const CharIntMap& symbol_table,
                                     const std::string& str);

  const std::vector<SymbolVector>& GetPatterns() const { return patterns_; }
  size_t PatternCount() const { return patterns_.at(0).size(); }
  size_t SequenceCount() const { return patterns_.size(); }
  size_t SiteCount() const { return alignment_.Length(); }
  const std::vector<double>& GetWeights() const { return weights_; }
  // Make a flattened partial likelihood vector for a given sequence, where anything
  // above 4 is given a uniform distribution.
  const std::vector<double> GetPartials(size_t sequence_idx) const;

  static SitePattern HelloSitePattern() {
    return SitePattern(Alignment::HelloAlignment(), {{PackInts(0, 1), "mars"},
                                                     {PackInts(1, 1), "saturn"},
                                                     {PackInts(2, 1), "jupiter"}});
  }

 private:
  Alignment alignment_;
  TagStringMap tag_taxon_map_;
  // The first index of patterns_ is across sequences, and the second is across site
  // patterns.
  std::vector<SymbolVector> patterns_;
  // The number of times each site pattern was seen in the alignment.
  std::vector<double> weights_;

  void Compress();
  static int SymbolTableAt(const CharIntMap& symbol_table, char c);
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("SitePattern") {
  CharIntMap symbol_table = SitePattern::GetSymbolTable();
  SymbolVector symbol_vector = SitePattern::SymbolVectorOf(symbol_table, "-tgcaTGCA?");
  SymbolVector correct_symbol_vector = {4, 3, 2, 1, 0, 3, 2, 1, 0, 4};
  CHECK_EQ(symbol_vector, correct_symbol_vector);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_SITE_PATTERN_HPP_
