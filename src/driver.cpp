// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "driver.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
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

// This parser will allow anything before the first '('.
TreeCollection::TreeCollectionPtr Driver::ParseNewick(std::ifstream &in) {
  yy::parser parser_instance(*this);

  parser_instance.set_debug_level(trace_parsing_);

  std::string line;
  unsigned int line_number = 1;
  Tree::TreePtrVector trees;
  while (std::getline(in, line)) {
    // Set the Bison location line number properly so we get useful error
    // messages.
    location_.initialize(nullptr, line_number);
    line_number++;
    auto tree_start = line.find_first_of("(");
    if (line.size() > 0 && tree_start != std::string::npos) {
      // Erase any characters before the first '('.
      line.erase(0, tree_start);
      trees.push_back(ParseString(&parser_instance, line));
    }
  }
  in.close();
  return std::make_shared<TreeCollection>(std::move(trees),
                                          this->TagTaxonMap());
}

TreeCollection::TreeCollectionPtr Driver::ParseNewickFile(
    const std::string &fname) {
  std::ifstream in(fname.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fname << std::endl;
    abort();
  }
  return ParseNewick(in);
}

TreeCollection::TreeCollectionPtr Driver::ParseNexusFile(
    const std::string &fname) {
  std::ifstream in(fname.c_str());
  try {
    if (!in) {
      throw std::runtime_error("Cannot open file.");
    }

    std::string line;
    std::getline(in, line);
    if (line != "#NEXUS") {
      throw std::runtime_error(
          "Putative Nexus file doesn't begin with #NEXUS.");
    }
    do {
      if (in.eof()) {
        throw std::runtime_error(
            "Finished reading and couldn't find 'begin trees;'");
        abort();
      }
      std::getline(in, line);
    } while (line != "begin trees;");
    std::getline(in, line);
    std::regex translate_start("^\\s*translate");
    if (!std::regex_match(line, translate_start)) {
      throw std::runtime_error("Missing translate block.");
    }
    std::getline(in, line);
    std::regex translate_item_regex("^\\s*(\\d+)\\s([^,]*)[,;]$");
    std::smatch match;
    while (std::regex_match(line, match, translate_item_regex)) {
      std::cout << match[1].str() << std::endl;
      std::cout << match[2].str() << std::endl;
      // Semicolon marks the end of the translate block.
      if (match[3].str() == ";") {
        break;
      }
      std::getline(in, line);
      if (in.eof()) {
        throw std::runtime_error(
            "Encountered EOF while parsing translate block.");
      }
    }
    std::cout << line << std::endl;
  } catch (const std::exception &exception) {
    std::cerr << "Problem parsing '" << fname << "':\n";
    std::cerr << exception.what() << std::endl;
    abort();
  }
  return ParseNewick(in);
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
