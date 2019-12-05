// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "site_pattern.hpp"
#include <cstdio>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "intpack.hpp"
#include "sugar.hpp"

// DNA assumption here.
CharIntMap SitePattern::GetSymbolTable() {
  CharIntMap table({{'A', 0},
                    {'C', 1},
                    {'G', 2},
                    {'T', 3},
                    {'a', 0},
                    {'c', 1},
                    {'g', 2},
                    {'t', 3},
                    {'-', 4},
                    {'N', 4},
                    {'X', 4},
                    {'?', 4},
                    // Treat degenerate nucleotides as gaps.
                    {'K', 4},
                    {'M', 4},
                    {'R', 4},
                    {'U', 4},
                    {'W', 4},
                    {'Y', 4}});
  return table;
}

int SitePattern::SymbolTableAt(const CharIntMap &symbol_table, char c) {
  auto search = symbol_table.find(c);
  if (search == symbol_table.end()) {
    char error[50];
    std::snprintf(error, sizeof(error), "Symbol '%c' not known.", c);
    Failwith(error);
  }
  return search->second;
}

SymbolVector SitePattern::SymbolVectorOf(const CharIntMap &symbol_table,
                                         const std::string &str) {
  SymbolVector v(str.size());
  for (size_t i = 0; i < str.size(); i++) {
    v[i] = SymbolTableAt(symbol_table, str[i]);
  }
  return v;
}

struct IntVectorHasher {
  int operator()(const std::vector<int> &values) const {
    int hash = values[0];
    for (size_t i = 1; i < values.size(); i++) {
      hash ^= values[i] + 0x9e3779b + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};

void SitePattern::Compress() {
  CharIntMap symbol_table = GetSymbolTable();
  size_t sequence_length = alignment_.Length();
  std::unordered_map<SymbolVector, double, IntVectorHasher> patterns;

  // Build an unordered map of patterns.
  for (size_t pos = 0; pos < sequence_length; pos++) {
    SymbolVector pattern(alignment_.SequenceCount());
    for (const auto &iter : tag_taxon_map_) {
      auto taxon_number = static_cast<size_t>(UnpackFirstInt(iter.first));
      auto symbol_to_find = alignment_.at(iter.second)[pos];
      pattern[taxon_number] = SymbolTableAt(symbol_table, symbol_to_find);
    }
    if (patterns.find(pattern) == patterns.end()) {
      SafeInsert(patterns, pattern, 1.);
    } else {
      patterns[pattern]++;
    }
  }

  // Collect the site patterns per taxon.
  for (const auto &iter_tag_taxon : tag_taxon_map_) {
    SymbolVector compressed_sequence;
    auto taxon_number = static_cast<size_t>(UnpackFirstInt(iter_tag_taxon.first));
    for (const auto &iter_patterns : patterns) {
      compressed_sequence.push_back(iter_patterns.first[taxon_number]);
    }
    patterns_[taxon_number] = compressed_sequence;
  }

  // Collect the site weights.
  for (const auto &iter : patterns) {
    weights_.push_back(iter.second);
  }
}
