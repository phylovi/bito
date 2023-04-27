#pragma once

#include "../src/alignment.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Alignment") {
  auto alignment = Alignment::ReadFasta("data/hello.fasta");
  CHECK_EQ(alignment, Alignment::HelloAlignment());
  CHECK(alignment.IsValid());

  CHECK_THROWS(alignment.ExtractSingleColumnAlignment(31));
  Alignment first_col_expected =
      Alignment({{"mars", "C"}, {"saturn", "G"}, {"jupiter", "G"}});
  CHECK_EQ(alignment.ExtractSingleColumnAlignment(0), first_col_expected);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
