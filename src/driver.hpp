// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_DRIVER_HPP_
#define SRC_DRIVER_HPP_
#include <map>
#include <string>
#include <vector>
#include "parser.hpp"
#include "tree_collection.hpp"
#include "typedefs.hpp"

// Give Flex the prototype of yylex we want ...
#define YY_DECL yy::parser::symbol_type yylex(Driver& drv)
// ... and declare it for the parser's sake.
YY_DECL;

// Conducting the scanning and parsing of trees.
// Note that this class is only for parsing a collection of trees on the same
// taxon set.
class Driver {
 public:
  Driver();

  // The next available id for parsing the first tree.
  uint32_t next_id_;
  // Is this the first tree we have parsed? The first tree gets to set up
  // indexing for the taxon names.
  bool first_tree_;
  // Whether to generate parser debug traces.
  bool trace_parsing_;
  // Whether to generate scanner debug traces.
  bool trace_scanning_;
  // The most recent tree parsed.
  Tree::TreePtr latest_tree_;
  // Map from taxon names to their numerical identifiers.
  std::map<std::string, uint32_t> taxa_;
  // The token's location, used by the scanner to give good debug info.
  TagDoubleMap branch_lengths_;
  // The token's location, used by the scanner to give good debug info.
  yy::location location_;

  // Clear out stored state.
  void Clear();
  // Scan a string with flex.
  void ScanString(const std::string& str);
  // Parse a string with an existing parser object.
  Tree::TreePtr ParseString(yy::parser* parser_instance,
                            const std::string& str);
  // Make a parser and then parse a string for a one-off parsing.
  TreeCollection::TreeCollectionPtr ParseString(const std::string& s);
  // Run the parser on a Newick stream.
  TreeCollection::TreeCollectionPtr ParseNewick(std::ifstream& in);
  // Run the parser on a Newick file.
  TreeCollection::TreeCollectionPtr ParseNewickFile(const std::string& fname);
  // Run the parser on a Nexus file.
  TreeCollection::TreeCollectionPtr ParseNexusFile(const std::string& fname);
  // Make the map from the edge tags of the tree to the taxon names from taxa_.
  TagStringMap TagTaxonMap();
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Driver") {
  Driver driver;

  std::vector<std::string> newicks = {
      "(a,b,c,d);",
      "((b,a),c);",
      "((a:1.1,b:2):0.4,c:3):0;",
      "(x,(a:1.1,(b:2,(quack:0.1,duck))),c:3):1.1;",
  };
  for (const auto& newick : newicks) {
    auto collection = driver.ParseString(newick);
    CHECK_EQ(newick, collection->Trees()[0]->Newick(collection->TagTaxonMap()));
  }
  driver.Clear();
  auto nexus_collection = driver.ParseNexusFile("data/DS1.subsampled.t");
  CHECK_EQ(nexus_collection->TreeCount(), 10);
  driver.Clear();
  auto newick_collection = driver.ParseNewickFile("data/DS1.subsampled.t.nwk");
  CHECK(nexus_collection == newick_collection);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_DRIVER_HPP_
