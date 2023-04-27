#pragma once

#include "../src/taxon_name_munging.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TaxonNameMunging") {
  using namespace TaxonNameMunging;

  std::string unquoted_test(R"raw(hello 'there" friend)raw");
  std::string double_quoted_test(R"raw("this is a \" test")raw");
  std::string double_quoted_dequoted(R"raw(this is a " test)raw");
  std::string single_quoted_test(R"raw('this is a \' test')raw");
  std::string single_quoted_dequoted(R"raw(this is a ' test)raw");

  CHECK_EQ(QuoteString(unquoted_test), R"raw("hello 'there\" friend")raw");
  CHECK_EQ(DequoteString(double_quoted_test), double_quoted_dequoted);
  CHECK_EQ(DequoteString(single_quoted_test), single_quoted_dequoted);
  CHECK_EQ(DequoteString(QuoteString(unquoted_test)), unquoted_test);

  TagStringMap test_map(
      {{2, unquoted_test}, {3, double_quoted_test}, {5, single_quoted_test}});
  TagStringMap expected_test_map(
      {{2, unquoted_test}, {3, double_quoted_dequoted}, {5, single_quoted_dequoted}});
  CHECK_EQ(expected_test_map, DequoteTagStringMap(test_map));

  // Test of TagDateMapOfTagTaxonMap appears in rooted_sbn_instance.hpp.
}
#endif  // DOCTEST_LIBRARY_INCLUDED
