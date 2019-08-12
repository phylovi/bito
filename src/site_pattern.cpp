// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "site_pattern.hpp"
#include <unordered_map>
#include <utility>
#include <vector>
#include "intpack.hpp"
#include "sugar.hpp"

// DNA assumption here.
CharIntMap GetSymbolTable() {
  CharIntMap table({{'A', 0},
                    {'C', 1},
                    {'G', 2},
                    {'T', 3},
                    {'a', 0},
                    {'c', 1},
                    {'g', 2},
                    {'t', 3},
                    {'-', 4}});
  return table;
}

struct VectorHasher {
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
  std::unordered_map<SymbolVector, double, VectorHasher> patterns;

  // Build an unordered map of patterns.
  for (size_t i = 0; i < sequence_length; i++) {
    SymbolVector pattern(alignment_.SequenceCount());
    for (const auto &iter : tag_taxon_map_) {
      size_t taxon_number = static_cast<size_t>(UnpackFirstInt(iter.first));
      pattern[taxon_number] = symbol_table.at(alignment_.at(iter.second)[i]);
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
    size_t taxon_number =
        static_cast<size_t>(UnpackFirstInt(iter_tag_taxon.first));
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
