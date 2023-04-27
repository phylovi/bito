#pragma once

#include "../src/mmapped_plv.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("MmappedNucleotidePLV") {
  MmappedNucleotidePLV mmapped_plv("_ignore/mmapped_plv.data", 10);
  auto plvs = mmapped_plv.Subdivide(2);
  for (const auto &plv : plvs) {
    CHECK_EQ(plv.rows(), MmappedNucleotidePLV::base_count_);
    CHECK_EQ(plv.cols(), 5);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
