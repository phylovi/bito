// Created by Mathieu Fourment on 23/7/19.
// Copyright © 2019 University of Technology Sydney. All rights reserved.
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
};

#endif  // SRC_SITE_PATTERN_HPP_
