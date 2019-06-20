#include "driver.hpp"

#include <fstream>
#include <iostream>

#include "parser.hpp"
#include "sbn.hpp"


Driver::Driver()
    : next_id_(0),
      first_tree_(true),
      trace_parsing_(false),
      trace_scanning_(false),
      latest_tree_(nullptr) {}


// TODO return the trees
void Driver::parse_file(const std::string &fname) {
  Node::NodePtr treePtr;

  yy::parser parserObject(*this);
  parserObject.set_debug_level(trace_parsing_);

  std::ifstream in(fname.c_str());
  if (!in) {
    // TODO do we want to raise an exception?
    std::cerr << "Cannot open the File : " << fname << std::endl;
    return;
  }
  std::string str;
  unsigned int line_number = 1;
  while (std::getline(in, str)) {
    location.initialize(nullptr, line_number);
    // Line contains string of length > 0 then save it in vector
    if (str.size() > 0) {
      std::cout << str << std::endl;
      treePtr = parse_string(parserObject, str);
      std::cout << treePtr->ToNewick() << std::endl;
      //  for (auto &x : taxa) {
      //    std::cout << x.first << " => " << x.second << '\n';
      //  }
      std::cout << std::endl;
    }
    line_number++;
  }
  in.close();
}


// Parse a string with an existing parser object.
Node::NodePtr Driver::parse_string(yy::parser &parserObject,
                                   const std::string &str) {
  this->scan_string(str);
  int return_code = parserObject();
  // TODO
  assert(return_code == 0);
  return latest_tree_;
}


// Make a parser and then parse a string for a one-off parsing.
Node::NodePtr Driver::parse_string(const std::string &str) {
  yy::parser parserObject(*this);
  parserObject.set_debug_level(trace_parsing_);
  return parse_string(parserObject, str);
}

// Note that a number of Driver methods are implemented in scanner.ll.
