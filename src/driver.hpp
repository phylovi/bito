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

  // TODO: reconsider result
  int result;
  // For indexing trees.
  int id_counter;
  // Is this the first tree we have parsed? The first tree gets to set up
  // indexing for the taxon names.
  bool first_tree;

  // Run the parser on file F.  Return 0 on success.
  int parse_file(const std::string& fname);
  // The name of the file being parsed.
  std::string fname;

  // Whether to generate parser debug traces.
  bool trace_parsing;
  // Whether to generate scanner debug traces.
  bool trace_scanning;

  int parse_string(const std::string& s);
  int parse_string(yy::parser& parserObject, const std::string& s);
  void scan_string(const std::string& s);

  // The token's location used by the scanner.
  yy::location location;
};
#endif  // ! DRIVER_HH
