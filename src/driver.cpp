#include "driver.hpp"

#include <fstream>
#include <iostream>
#include <memory>

#include "parser.hpp"
#include "sbn.hpp"


Driver::Driver()
    : next_id_(0),
      first_tree_(true),
      trace_parsing_(false),
      trace_scanning_(false),
      latest_tree_(nullptr) {}


Node::NodePtrVecPtr Driver::ParseFile(const std::string &fname) {
  Node::NodePtr tree;
  yy::parser parser_instance(*this);

  parser_instance.set_debug_level(trace_parsing_);

  std::ifstream in(fname.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fname << std::endl;
    abort();
  }
  std::string line;
  unsigned int line_number = 1;
  auto trees = std::make_shared<Node::NodePtrVec>();
  while (std::getline(in, line)) {
    // Set the Bison location line number properly so we get useful error
    // messages.
    location_.initialize(nullptr, line_number);
    line_number++;
    if (line.size() > 0) {
      tree = ParseString(&parser_instance, line);
      trees->push_back(tree);
    }
  }
  in.close();
  return trees;
}

Node::NodePtr Driver::ParseString(yy::parser *parser_instance,
                                  const std::string &str) {
  this->ScanString(str);
  int return_code = (*parser_instance)();
  if (return_code != 0) {
    abort();
  }
  return latest_tree_;
}

Node::NodePtr Driver::ParseString(const std::string &str) {
  yy::parser parser_instance(*this);
  parser_instance.set_debug_level(trace_parsing_);
  return ParseString(&parser_instance, str);
}

// Note that a number of Driver methods are implemented in scanner.ll.
