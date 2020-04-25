// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// This class drives tree parsing.

#ifndef SRC_DRIVER_HPP_
#define SRC_DRIVER_HPP_
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "parser.hpp"
#include "rooted_tree_collection.hpp"
#include "sugar.hpp"
#include "tree_collection.hpp"
#include "unrooted_tree_collection.hpp"

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

  // These member variables are public so that the parser and lexer can access them.
  // Nevertheless, the reasonable thing to do is to interact with this class from the
  // outside using the public method interface below.

  // The next available id for parsing the first tree.
  uint32_t next_id_;
  // Do we already have the taxon names in taxa_? If not, they get initialized with the
  // first tree parsed.
  bool taxa_complete_;
  // Debug level for parser.
  int trace_parsing_;
  // Whether to generate scanner debug traces.
  bool trace_scanning_;
  // The most recent tree parsed.
  std::shared_ptr<Tree> latest_tree_;
  // Map from taxon names to their numerical identifiers.
  std::map<std::string, uint32_t> taxa_;
  // The token's location, used by the scanner to give good debug info.
  TagDoubleMap branch_lengths_;
  // The token's location, used by the scanner to give good debug info.
  yy::location location_;

  // These three parsing methods also remove quotes from Newick strings and Nexus files.
  // Make a parser and then parse a string for a one-off parsing.
  TreeCollection ParseString(const std::string& s);
  // Run the parser on a Newick file.
  TreeCollection ParseNewickFile(const std::string& fname);
  // Run the parser on a Nexus file.
  TreeCollection ParseNexusFile(const std::string& fname);
  // Clear out stored state.
  void Clear();
  // Make the map from the edge tags of the tree to the taxon names from taxa_.
  TagStringMap TagTaxonMap();

 private:
  // Scan a string with flex.
  void ScanString(const std::string& str);
  // Parse a string with an existing parser object.
  Tree ParseString(yy::parser* parser_instance, const std::string& str);
  // Run the parser on a Newick stream.
  TreeCollection ParseNewick(std::ifstream& in);
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Driver") {
  Driver driver;

  std::vector<std::string> newicks = {
      "(a:0,b:0,c:0,d:0):0;",
      "((b:0,a:0):0,c:0):0;",
      "((a:1.1,b:2):0.4,c:3):0;",
      "(x:0,(a:1.1,(b:2,(quack:0.1,duck:0):0):0):0,c:3):1.1;",
  };
  for (const auto& newick : newicks) {
    auto collection = driver.ParseString(newick);
    CHECK_EQ(newick, collection.Trees()[0].Newick(collection.TagTaxonMap()));
  }
  driver.Clear();
  // Note that the order of the taxa is given by the order in the translate table, not
  // by the short names. We use that here to make sure that the ordering of the taxa is
  // the same as that in the newick file below so that they can be compared.
  auto nexus_collection = driver.ParseNexusFile("data/DS1.subsampled_10.t.reordered");
  CHECK_EQ(nexus_collection.TreeCount(), 10);
  driver.Clear();
  auto newick_collection = driver.ParseNewickFile("data/DS1.subsampled_10.t.nwk");
  CHECK_EQ(nexus_collection, newick_collection);
  driver.Clear();
  auto five_taxon = driver.ParseNewickFile("data/five_taxon.nwk");
  std::vector<std::string> correct_five_taxon_names({"x0", "x1", "x2", "x3", "x4"});
  CHECK_EQ(five_taxon.TaxonNames(), correct_five_taxon_names);
  // Check that we can parse BEAST trees with [&comments], and that the different
  // formatting of the translate block doesn't trip us up.
  auto beast_nexus = driver.ParseNexusFile("data/test_beast_tree_parsing.nexus");
  // These are the taxa, in order, taken directly from the nexus file:
  StringVector beast_taxa = {
      "aDuckA_1976",     "aDuckB_1977",   "aItaly_1987",   "aMallard_1985",
      "hCHR_1983",       "hCambr_1939",   "hFortMon_1947", "hKiev_1979",
      "hLenin_1954",     "hMongol_1985",  "hMongol_1991",  "hNWS_1933",
      "hPR_1934",        "hSCar_1918.00", "hScot_1994",    "hSuita_1989",
      "hUSSR_1977",      "sEhime_1980",   "sIllino_1963",  "sIowa_1930",
      "sNebrask_1992",   "sNewJers_1976", "sStHya_1991",   "sWiscons_1961",
      "sWiscons_1.998e3"};
  CHECK_EQ(beast_nexus.TaxonNames(), beast_taxa);
  // Check that we got the whole tree.
  for (const auto& [topology, count] : beast_nexus.TopologyCounter()) {
    CHECK_EQ(topology->LeafCount(), beast_taxa.size());
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_DRIVER_HPP_
