// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "driver.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <unordered_map>
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
    auto previous_position = in.tellg();
    std::unordered_map<std::string, std::string> translator;
    while (std::regex_match(line, match, translate_item_regex)) {
      assert(translator.insert({match[1].str(), match[2].str()}).second);
      // Semicolon marks the end of the translate block.
      if (match[3].str() == ";") {
        break;
      }
      previous_position = in.tellg();
      std::getline(in, line);
      if (in.eof()) {
        throw std::runtime_error(
            "Encountered EOF while parsing translate block.");
      }
    }
    // Back up one line to hit the first tree.
    in.seekg(previous_position);
    // Now we make a new TagTaxonMap to replace the one with numbers in place of
    // taxon names.
    auto pre_translation = ParseNewick(in);
    TagStringMap translated_taxon_map;
    for (const auto &iter : pre_translation->TagTaxonMap()) {
      auto search = translator.find(iter.second);
      if (search == translator.end()) {
        throw std::runtime_error("Couldn't find a translation table entry.");
      }
      assert(translated_taxon_map.insert({iter.first, search->second}).second);
    }
    return std::make_shared<TreeCollection>(std::move(pre_translation->Trees()),
                                            std::move(translated_taxon_map));
  } catch (const std::exception &exception) {
    std::cerr << "\nProblem parsing '" << fname << "':\n";
    std::cerr << exception.what() << std::endl;
    abort();
  }
}

Tree::TreePtr Driver::ParseString(yy::parser *parser_instance,
                                  const std::string &str) {
  // Scan the string using the lexer into hidden state.
  this->ScanString(str);
  // Parse the scanned string.
  int return_code = (*parser_instance)();
  if (return_code != 0) {
    std::cout << "Parser had nonzero return value.\n";
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
