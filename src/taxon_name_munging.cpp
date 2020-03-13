#include "taxon_name_munging.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

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
