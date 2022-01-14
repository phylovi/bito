// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "csv.h"
#include "sugar.hpp"

namespace CSV {

StringDoubleMap StringDoubleMapOfCSV(const std::string& csv_path);
void StringDoubleVectorToCSV(const StringDoubleVector& v, const std::string& csv_path);

}  // namespace CSV

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("CSV I/O") {
  std::string csv_test_file_path = "_ignore/for_csv_test.csv";
  StringDoubleVector input = {{"hi", 1e9}, {"lo", -4.}};
  CSV::StringDoubleVectorToCSV(input, csv_test_file_path);
  auto result = CSV::StringDoubleMapOfCSV(csv_test_file_path);
  StringDoubleMap correct_result = {{"hi", 1e9}, {"lo", -4.}};
  CHECK_EQ(result, correct_result);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
