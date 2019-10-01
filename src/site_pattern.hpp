// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SITE_PATTERN_HPP_
#define SRC_SITE_PATTERN_HPP_

#include <vector>
#include "alignment.hpp"
#include "sugar.hpp"

class SitePattern {
 public:
  SitePattern() {}
  explicit SitePattern(Alignment& alignment, const TagStringMap& tag_taxon_map)
      : alignment_(alignment), tag_taxon_map_(tag_taxon_map) {
    patterns_.resize(alignment.SequenceCount());
    Compress();
  }

  static CharIntMap GetSymbolTable();
  static SymbolVector SymbolVectorOf(const std::string& str,
                                     const CharIntMap& symbol_table);

  const std::vector<SymbolVector>& GetPatterns() const { return patterns_; }
  size_t PatternCount() const { return patterns_.at(0).size(); }
  size_t SequenceCount() const { return patterns_.size(); }
  const std::vector<double>& GetWeights() const { return weights_; }

 private:
  Alignment alignment_;
  TagStringMap tag_taxon_map_;
  std::vector<SymbolVector> patterns_;
  std::vector<double> weights_;

  void Compress();
  static void FailwithUnknownSymbol(char c);
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("SitePattern") {
  CharIntMap symbol_table = SitePattern::GetSymbolTable();
  SymbolVector symbol_vector =
      SitePattern::SymbolVectorOf("-tgcaTGCA?", symbol_table);
  SymbolVector correct_symbol_vector = {4, 3, 2, 1, 0, 3, 2, 1, 0, 4};
  CHECK_EQ(symbol_vector, correct_symbol_vector);

}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_SITE_PATTERN_HPP_
