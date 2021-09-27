// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <string>

#include "rooted_sbn_instance.hpp"
#include "stick_breaking_transform.hpp"
#include "taxon_name_munging.hpp"
#include "unrooted_sbn_instance.hpp"

// NOTE: This file is automatically generated from `test/prep/doctest.py`. Don't edit!

TEST_CASE("Node") {
  Driver driver;

  std::vector<std::string> trace;
  auto t =
      driver.ParseString("((((0_1,1_1),(2_1,3_1)),4_1),((5_1,(6_1,7_1)),(8_1,9_1)));")
          .Trees()[0];

  // preorder:
  t.Topology()->Preorder(
      [&trace](const Node* node) { trace.push_back(node->TagString()); });
  CHECK(std::vector<std::string>({"9_10", "4_5", "3_4", "1_2", "0_1", "1_1", "3_2",
                                  "2_1", "3_1", "4_1", "9_5", "7_3", "5_1", "7_2",
                                  "6_1", "7_1", "9_2", "8_1", "9_1"}) == trace);
  trace.clear();

  // postorder:
  t.Topology()->Postorder(
      [&trace](const Node* node) { trace.push_back(node->TagString()); });
  CHECK(std::vector<std::string>({"0_1", "1_1", "1_2", "2_1", "3_1", "3_2", "3_4",
                                  "4_1", "4_5", "5_1", "6_1", "7_1", "7_2", "7_3",
                                  "8_1", "9_1", "9_2", "9_5", "9_10"}) == trace);
  trace.clear();

  // levelorder:
  t.Topology()->LevelOrder(
      [&trace](const Node* node) { trace.push_back(node->TagString()); });
  CHECK(std::vector<std::string>({"9_10", "4_5", "9_5", "3_4", "4_1", "7_3", "9_2",
                                  "1_2", "3_2", "5_1", "7_2", "8_1", "9_1", "0_1",
                                  "1_1", "2_1", "3_1", "6_1", "7_1"}) == trace);
  trace.clear();
}
