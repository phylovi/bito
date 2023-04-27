// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This class drives tree parsing.

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "parser.hpp"
#include "sugar.hpp"
#include "tree_collection.hpp"

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
  // Run the parser on a gzip-ed Newick file.
  TreeCollection ParseNewickFileGZ(const std::string& fname);
  // Run the parser on a Nexus file. The Nexus file must have a translate block, and the
  // leaf tags are assigned according to the order of names in the translate block.
  TreeCollection ParseNexusFile(const std::string& fname);
  // Run the parser on a gzip-ed Nexus file. Check ParseNexusFile() for details.
  TreeCollection ParseNexusFileGZ(const std::string& fname);
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
  TreeCollection ParseNewick(std::istream& in);
  // Runs ParseNewick() and dequotes the resulting trees.
  TreeCollection ParseAndDequoteNewick(std::istream& in);
  // Run the parser on a Nexus stream.
  TreeCollection ParseNexus(std::istream& in);
};
