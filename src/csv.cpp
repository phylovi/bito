// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "csv.hpp"

#include <fstream>

// Given a headerless 2-column CSV of quoted string keys then double values, this parses
// the CSV into a StringDoubleMap.
StringDoubleMap CSV::StringDoubleMapOfCSV(const std::string& in_path) {
  io::CSVReader<2, io::trim_chars<' ', '\t'>, io::double_quote_escape<',', '"'>> csv_in(
      in_path);
  std::string key;
  double value;
  StringDoubleMap string_double_map;
  while (csv_in.read_row(key, value)) {
    auto search = string_double_map.find(key);
    if (search == string_double_map.end()) {
      string_double_map.insert({key, value});
    } else {
      Failwith("Key " + key + " found twice in " + in_path +  // NOLINT
               "when turning it into a map.");
    }
  }
  return string_double_map;
}

void CSV::StringDoubleVectorToCSV(const StringDoubleVector& vect,
                                  const std::string& out_path) {
  std::ofstream out_stream(out_path);
  for (const auto& [s, value] : vect) {
    out_stream << s << "," << value << std::endl;
  }
  if (out_stream.bad()) {
    Failwith("Failure writing to " + out_path);
  }
  out_stream.close();
}
