// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "driver.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <utility>

#include "parser.hpp"
#include "tree.hpp"


Driver::Driver()
    : next_id_(0),
      first_tree_(true),
      trace_parsing_(false),
      trace_scanning_(false),
      latest_tree_(nullptr) {}

void Driver::Clear() {
  next_id_ = 0;
  first_tree_ = true;
  trace_parsing_ = false;
  trace_scanning_ = false;
  latest_tree_ = nullptr;
  branch_lengths_.clear();
}

Tree::TreePtrCounterPtr Driver::ParseFile(const std::string &fname) {
  Tree::TreePtr tree;
  yy::parser parser_instance(*this);

  parser_instance.set_debug_level(trace_parsing_);

  std::ifstream in(fname.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fname << std::endl;
    abort();
  }
  std::string line;
  unsigned int line_number = 1;
  auto trees = std::make_shared<Tree::TreePtrCounter>();
  while (std::getline(in, line)) {
    // Set the Bison location line number properly so we get useful error
    // messages.
    location_.initialize(nullptr, line_number);
    line_number++;
    if (line.size() > 0) {
      tree = ParseString(&parser_instance, line);
      auto search = trees->find(tree);
      if (search == trees->end()) {
        assert(trees->insert(std::make_pair(tree, 1)).second);
      } else {
        search->second++;
      }
    }
  }
  in.close();
  return trees;
}

Tree::TreePtr Driver::ParseString(yy::parser *parser_instance,
                                  const std::string &str) {
  this->ScanString(str);
  int return_code = (*parser_instance)();
  if (return_code != 0) {
    abort();
  }
  return latest_tree_;
}

Tree::TreePtr Driver::ParseString(const std::string &str) {
  yy::parser parser_instance(*this);
  parser_instance.set_debug_level(trace_parsing_);
  return ParseString(&parser_instance, str);
}

// Note that a number of Driver methods are implemented in scanner.ll.
