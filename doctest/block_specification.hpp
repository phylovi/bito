#pragma once

#include "../src/block_specification.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("BlockSpecification") {
  // As an example, kazoo has 4 parameters, and jordan has 23.
  BlockSpecification spec({{"kazoo", 4}, {"jordan", 23}});
  // The specification stores the starting index and then the number of
  // parameters. Because we're using an ordered map, jordan has a lower index
  // than kazoo.
  const auto correct_spec_map = BlockSpecification::UnderlyingMapType(
      {{"entire", {0, 27}}, {"jordan", {0, 23}}, {"kazoo", {23, 4}}});
  CHECK_EQ(spec.GetMap(), correct_spec_map);
  spec.Append("entire turbo and boost",
              BlockSpecification({{"boost", 42}, {"turbo", 666}}));
  // Then after appending, the new stuff gets shifted down. For example, we find
  // boost at 23+4=27 and turbo at 27+42=69.
  auto correct_appended_map = BlockSpecification::UnderlyingMapType(
      {{"boost", {27, 42}},                    // 23+4=27
       {"entire", {0, 735}},                   // 4+23+42+666=735
       {"entire turbo and boost", {27, 708}},  // 42+666=708
       {"jordan", {0, 23}},                    //
       {"kazoo", {23, 4}},                     //
       {"turbo", {69, 666}}});                 // 27+42=69
  CHECK_EQ(spec.GetMap(), correct_appended_map);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
