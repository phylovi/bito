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

TagDoubleMap TaxonNameMunging::TagDateMapOfTagTaxonMap(TagStringMap tag_taxon_map) {
  TagDoubleMap tag_date_map;
  std::regex date_regex(R"raw(^.+_(\d*\.?\d+(?:[eE][-+]?\d+)?)$)raw");
  std::smatch match_date;
  bool first_pass_through_parsing_loop = true;
  bool have_parsed_a_date = false;
  double max_date = DOUBLE_NEG_INF;
  for (auto &[tag, taxon] : tag_taxon_map) {
    if (std::regex_match(taxon, match_date, date_regex)) {
      double date = std::stod(match_date[1].str());
      max_date = std::max(date, max_date);
      SafeInsert(tag_date_map, tag, date);
      if (first_pass_through_parsing_loop) {
        have_parsed_a_date = true;
      } else {
        Assert(have_parsed_a_date,
               "We couldn't parse dates for a while, but we could parse:" + taxon);
      }
    } else {  // We couldn't parse a date.
      if (!first_pass_through_parsing_loop && have_parsed_a_date) {
        Failwith("We did parse at least one date, but couldn't parse:" + taxon);
      }
      // This is the first pass through the loop or we haven't parsed a date previously.
      // In this case we fill tag_date_map with zeroes.
      SafeInsert(tag_date_map, tag, 0.);
    }
    first_pass_through_parsing_loop = false;
  }
  if (have_parsed_a_date) {
    for (auto &[id, date] : tag_date_map) {
      date = max_date - date;
    }
  }
  return tag_date_map;
}
