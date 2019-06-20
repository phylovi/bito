#include "driver.hpp"
#include "parser.hpp"

#include <fstream>
#include <iostream>

driver::driver()
    : id_counter(0),
      first_tree(false),
      trace_parsing(false),
      trace_scanning(false) {}


int driver::parse_file(const std::string &f) {
  int return_code;

  fname = f;
  yy::parser parserObject(*this);
  parserObject.set_debug_level(trace_parsing);

  std::ifstream in(fname.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fname << std::endl;
    return false;
  }
  std::string str;
  unsigned int line_number = 1;
  while (std::getline(in, str)) {
    location.initialize(nullptr, line_number);
    // Line contains string of length > 0 then save it in vector
    if (str.size() > 0) {
      return_code = parse_string(parserObject, str);
      if (return_code != 0) {
        std::cout << "problem parsing!\n";
        break;
      }
      std::cout << result << '\n';
      for (auto &x : taxa) {
        std::cout << x.first << " => " << x.second << '\n';
      }
    }
    line_number++;
  }
  in.close();
  return return_code;
}


// Parse a string with an existing parser object.
int driver::parse_string(yy::parser &parserObject, const std::string &str) {
  this->scan_string(str);
  int res = parserObject();
  return res;
}


// Make a parser and then parse a string for a one-off parsing.
int driver::parse_string(const std::string &str) {
  yy::parser parserObject(*this);
  parserObject.set_debug_level(trace_parsing);
  return parse_string(parserObject, str);
}

// Note that a number of driver methods are implemented in scanner.ll.
