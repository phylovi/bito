// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <functional>

#include "sugar.hpp"

namespace TaxonNameMunging {
std::string QuoteString(const std::string &in_str);
std::string DequoteString(const std::string &in_str);
TagStringMap DequoteTagStringMap(const TagStringMap &tag_string_map);

// Make each date in tag_date_map the given date minus the maximum date.
void MakeDatesRelativeToMaximum(TagDoubleMap &tag_date_map);

// Returns a map from each tag to zero.
TagDoubleMap ConstantDatesForTagTaxonMap(TagStringMap tag_taxon_map);

// Parses dates as numbers appearing after an underscore. Returns a map specifying the
// height of each taxon compared to the maximum date.
TagDoubleMap ParseDatesFromTagTaxonMap(TagStringMap tag_taxon_map);
}  // namespace TaxonNameMunging
