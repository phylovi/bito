#ifndef DRIVER_HH
#define DRIVER_HH
#include <map>
#include <string>
#include "parser.hpp"

// Give Flex the prototype of yylex we want ...
#define YY_DECL yy::parser::symbol_type yylex(driver& drv)
// ... and declare it for the parser's sake.
YY_DECL;

// Conducting the scanning and parsing of trees.
// Note that this class is only for parsing a collection of trees on the same
// taxon set.
class driver {
 public:
  driver();

  std::map<std::string, int> taxa;

  // The next available id for tree parsing.
  int next_id;
  // Is this the first tree we have parsed? The first tree gets to set up
  // indexing for the taxon names.
  bool first_tree;
  // Whether to generate parser debug traces.
  bool trace_parsing;
  // Whether to generate scanner debug traces.
  bool trace_scanning;
  // The most recent tree parsed.
  Node::NodePtr latest_tree;

  void scan_string(const std::string& s);
  Node::NodePtr parse_string(const std::string& s);
  Node::NodePtr parse_string(yy::parser& parserObject, const std::string& s);

  // Run the parser on file F.  Return 0 on success.
  void parse_file(const std::string& fname);

  // The token's location used by the scanner.
  yy::location location;
};
#endif  // ! DRIVER_HH
