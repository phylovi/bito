#ifndef DRIVER_HH
#define DRIVER_HH
#include <map>
#include <string>
#include "parser.hpp"

// Give Flex the prototype of yylex we want ...
#define YY_DECL yy::parser::symbol_type yylex(driver& drv)
// ... and declare it for the parser's sake.
YY_DECL;

// Conducting the whole scanning and parsing of Calc++.
class driver {
 public:
  driver();

  std::map<std::string, int> taxa;

  int result;
  int id_counter;

  void clear();

  // Run the parser on file F.  Return 0 on success.
  int parse_file(const std::string& fname);
  // The name of the file being parsed.
  std::string fname;

  // Whether to generate parser debug traces.
  bool trace_parsing;
  // Whether to generate scanner debug traces.
  bool trace_scanning;

  int parse_string(const std::string& s);
  void scan_string(const std::string& s);

  // The token's location used by the scanner.
  yy::location location;
};
#endif  // ! DRIVER_HH
