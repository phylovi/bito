#include "driver.hh"
#include "parser.hh"

#include <fstream>
#include <iostream>

driver::driver() : trace_parsing(false), trace_scanning(false) {
  id_counter = 0;
}

void driver::clear() {
  id_counter = 0;
  taxa.clear();
}

int driver::parse_file(const std::string &f) {
  int return_code;

  fname = f;
  location.initialize(&fname);
  yy::parser parserObject(*this);
  parserObject.set_debug_level(trace_parsing);

  std::ifstream in(fname.c_str());
  if(!in)
  {
  std::cerr << "Cannot open the File : "<<fname<<std::endl;
  return false;
  }
  std::string str;
  while (std::getline(in, str))
  {
  // Line contains string of length > 0 then save it in vector
  if(str.size() > 0) {
    this->scan_string(str);
    return_code = parserObject();
    if(return_code != 0) {
      std::cout << "problem parsing!\n";
      break;
      }
    std::cout << result << '\n';
    for (auto& x: taxa) {
      std::cout << x.first << " => " << x.second << '\n';
    }
    clear();
    }
  }
  in.close();
  return return_code;
}

int driver::parse_string(const std::string &s) {
  yy::parser parserObject(*this);
  parserObject.set_debug_level(trace_parsing);
  this->scan_string(s);
  int res = parserObject();
  return res;
}

// Note that a number of driver methods are implemented in scanner.ll.
