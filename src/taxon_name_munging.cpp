// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "taxon_name_munging.hpp"

#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>

#include "numerical_utils.hpp"

std::string TaxonNameMunging::QuoteString(const std::string &in_str) {
  std::stringstream ss;
  ss << std::quoted(in_str);
  return ss.str();
}

std::string TaxonNameMunging::DequoteString(const std::string &in_str) {
  if (in_str.empty()) {
    return std::string();
  }
  char delimiter = in_str.at(0);
  if (delimiter != '\'' && delimiter != '"') {
    return std::string(in_str);
  }  // else
  std::stringstream ss(in_str);
  std::string out_str;
  ss >> std::quoted(out_str, delimiter);
  return out_str;
}

TagStringMap TransformStringValues(std::function<std::string(const std::string &)> f,
                                   const TagStringMap &in_map) {
  TagStringMap out_map;
  for (const auto &[tag, value] : in_map) {
    out_map.insert({tag, f(value)});
  }
  return out_map;
}

TagStringMap TaxonNameMunging::DequoteTagStringMap(const TagStringMap &tag_string_map) {
  return TransformStringValues(DequoteString, tag_string_map);
}

void TaxonNameMunging::MakeDatesRelativeToMaximum(TagDoubleMap &tag_date_map) {
  double max_date = DOUBLE_NEG_INF;
  for (const auto &[_, date] : tag_date_map) {
    std::ignore = _;
    max_date = std::max(date, max_date);
  }
  for (auto &[id, date] : tag_date_map) {
    std::ignore = id;
    date = max_date - date;
  }
}

TagDoubleMap TaxonNameMunging::ConstantDatesForTagTaxonMap(TagStringMap tag_taxon_map) {
  TagDoubleMap tag_date_map;
  for (auto &[tag, taxon] : tag_taxon_map) {
    std::ignore = taxon;
    SafeInsert(tag_date_map, tag, 0.);
  }
  return tag_date_map;
}

TagDoubleMap TaxonNameMunging::ParseDatesFromTagTaxonMap(TagStringMap tag_taxon_map) {
  TagDoubleMap tag_date_map;
  std::regex date_regex(R"raw(^.+_(\d*\.?\d+(?:[eE][-+]?\d+)?)$)raw");
  std::smatch match_date;
  for (auto &[tag, taxon] : tag_taxon_map) {
    if (std::regex_match(taxon, match_date, date_regex)) {
      double date = std::stod(match_date[1].str());
      SafeInsert(tag_date_map, tag, date);
    } else {
      Failwith("Couldn't parse a date from:" + taxon);
    }
  }
  MakeDatesRelativeToMaximum(tag_date_map);
  return tag_date_map;
}
