// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "driver.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>
#include "parser.hpp"

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
  taxa_.clear();
  branch_lengths_.clear();
}

TreeCollection::TreeCollectionPtr Driver::ParseFile(const std::string &fname) {
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
  Tree::TreePtrVector trees;
  while (std::getline(in, line)) {
    // Set the Bison location line number properly so we get useful error
    // messages.
    location_.initialize(nullptr, line_number);
    line_number++;
    if (line.size() > 0) {
      trees.push_back(ParseString(&parser_instance, line));
    }
  }
  in.close();
  return std::make_shared<TreeCollection>(std::move(trees),
                                          this->TagTaxonMap());
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

TreeCollection::TreeCollectionPtr Driver::ParseString(const std::string &str) {
  Clear();
  yy::parser parser_instance(*this);
  parser_instance.set_debug_level(trace_parsing_);
  Tree::TreePtrVector trees = {ParseString(&parser_instance, str)};
  return std::make_shared<TreeCollection>(trees, this->TagTaxonMap());
}

TagStringMap Driver::TagTaxonMap() {
  TagStringMap m;
  for (const auto &iter : taxa_) {
    // These are leaves, so the number of leaves below is 1.
    m[PackInts(iter.second, 1)] = iter.first;
  }
  return m;
}

// Note that a number of Driver methods are implemented in scanner.ll.
