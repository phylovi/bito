#ifndef DRIVER_HH
#define DRIVER_HH
#include <map>
#include <string>
#include "parser.hh"

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
  int parse(const std::string& f);
  // The name of the file being parsed.
  std::string file;
  // Whether to generate parser debug traces.
  bool trace_parsing;
  int parse_string(const std::string& s);
  void scan_string(const std::string &s);

  // Handling the scanner.
  void scan_begin();
  void scan_end();
  // Whether to generate scanner debug traces.
  bool trace_scanning;
  // The token's location used by the scanner.
  yy::location location;
};
#endif  // ! DRIVER_HH
