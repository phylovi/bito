// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_DRIVER_HPP_
#define SRC_DRIVER_HPP_
#include <map>
#include <string>
#include "parser.hpp"

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
  int next_id_;
  // Is this the first tree we have parsed? The first tree gets to set up
  // indexing for the taxon names.
  bool first_tree_;
  // Whether to generate parser debug traces.
  bool trace_parsing_;
  // Whether to generate scanner debug traces.
  bool trace_scanning_;
  // The most recent tree parsed.
  Node::NodePtr latest_tree_;
  // Map from taxon names to their numerical identifiers.
  std::map<std::string, int> taxa_;
  // The token's location, used by the scanner to give good debug info.
  Tree::BranchLengthMap branch_lengths_;
  // The token's location, used by the scanner to give good debug info.
  yy::location location_;

  // Clear out stored state.
  void Clear();
  // Scan a string with flex.
  void ScanString(const std::string& str);
  // Parse a string with an existing parser object.
  Node::NodePtr ParseString(yy::parser* parser_instance,
                            const std::string& str);
  // Make a parser and then parse a string for a one-off parsing.
  Node::NodePtr ParseString(const std::string& s);
  // Run the parser on a file.
  Node::NodePtrCounterPtr ParseFile(const std::string& fname);
};
#endif  // SRC_DRIVER_HPP_
