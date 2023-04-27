// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "csv.h"
#include "sugar.hpp"

namespace CSV {

StringDoubleMap StringDoubleMapOfCSV(const std::string& csv_path);
void StringDoubleVectorToCSV(const StringDoubleVector& v, const std::string& csv_path);

}  // namespace CSV
