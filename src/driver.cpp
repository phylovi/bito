// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

// Based on
// https://www.gnu.org/software/bison/manual/html_node/Calc_002b_002b-Parsing-Driver.html#Calc_002b_002b-Parsing-Driver

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
      taxa_complete_(false),
      trace_parsing_(0),
      trace_scanning_(false),
      latest_tree_(nullptr) {}

void Driver::Clear() {
  next_id_ = 0;
  taxa_complete_ = false;
  trace_parsing_ = 0;
  trace_scanning_ = false;
  latest_tree_ = nullptr;
  taxa_.clear();
  branch_lengths_.clear();
}

// This parser will allow anything before the first '('.
TreeCollection Driver::ParseNewick(std::ifstream &in) {
  yy::parser parser_instance(*this);
  parser_instance.set_debug_level(trace_parsing_);
  std::string line;
  unsigned int line_number = 1;
  Tree::TreeVector trees;
  while (std::getline(in, line)) {
    // Set the Bison location line number properly so we get useful error
    // messages.
    location_.initialize(nullptr, line_number);
    line_number++;
    auto tree_start = line.find_first_of('(');
    if (!line.empty() && tree_start != std::string::npos) {
      // Erase any characters before the first '('.
      line.erase(0, tree_start);
      trees.push_back(ParseString(&parser_instance, line));
    }
  }
  in.close();
  return TreeCollection(std::move(trees), this->TagTaxonMap());
}

TreeCollection Driver::ParseNewickFile(const std::string &fname) {
  Clear();
  std::ifstream in(fname.c_str());
  if (!in) {
    Failwith("Cannot open the File : " + fname);
  }
  return ParseNewick(in);
}

TreeCollection Driver::ParseNexusFile(const std::string &fname) {
  Clear();
  std::ifstream in(fname.c_str());
  try {
    if (!in) {
      throw std::runtime_error("Cannot open file.");
    }
    std::string line;
    std::getline(in, line);
    if (line != "#NEXUS") {
      throw std::runtime_error("Putative Nexus file doesn't begin with #NEXUS.");
    }
    do {
      if (in.eof()) {
        throw std::runtime_error("Finished reading and couldn't find 'begin trees;'");
      }
      std::getline(in, line);
      // BEAST uses "Begin trees;" so we tolower here.
      line[0] = std::tolower(line[0]);
    } while (line != "begin trees;");
    std::getline(in, line);
    std::regex translate_start("^\\s*[Tt]ranslate");
    if (!std::regex_match(line, translate_start)) {
      throw std::runtime_error("Missing translate block.");
    }
    std::getline(in, line);
    std::regex translate_item_regex(R"raw(^\s*(\d+)\s([^,;]*)[,;]?$)raw");
    std::regex lone_semicolon_regex(R"raw(\s*;$)raw");
    std::smatch match;
    auto previous_position = in.tellg();
    TagStringMap long_name_taxon_map;
    uint32_t leaf_id = 0;
    while (std::regex_match(line, match, translate_item_regex)) {
      const auto short_name = match[1].str();
      const auto long_name = match[2].str();
      // We prepare taxa_ so that it can parse the short taxon names.
      SafeInsert(taxa_, short_name, leaf_id);
      // However, we keep the long names for the TagTaxonMap.
      SafeInsert(long_name_taxon_map, PackInts(leaf_id, 1), long_name);
      leaf_id++;
      // Semicolon marks the end of the translate block.
      // It appears at the end of a translation statement line in MrBayes.
      if (match[3].str() == ";") {
        break;
      }
      previous_position = in.tellg();
      std::getline(in, line);
      // BEAST has the ending semicolon on a line of its own.
      if (std::regex_match(line, match, lone_semicolon_regex)) {
        break;
      }
      if (in.eof()) {
        throw std::runtime_error("Encountered EOF while parsing translate block.");
      }
    }
    Assert(leaf_id > 0, "No taxa found in translate block!");
    taxa_complete_ = true;
    // Back up one line to hit the first tree.
    in.seekg(previous_position);
    // Now we make a new TagTaxonMap to replace the one with numbers in place of
    // taxon names.
    auto short_name_tree_collection = ParseNewick(in);
    // We're using the public member directly rather than the const accessor because we
    // want to move.
    return TreeCollection(std::move(short_name_tree_collection.trees_),
                          std::move(long_name_taxon_map));
  } catch (const std::exception &exception) {
    Failwith("Problem parsing '" + fname + "':\n" + exception.what());
  }
}

Tree Driver::ParseString(yy::parser *parser_instance, const std::string &str) {
  // Scan the string using the lexer into hidden state.
  this->ScanString(str);
  // Parse the scanned string.
  int return_code = (*parser_instance)();
  Assert(return_code == 0, "Parser had nonzero return value.");
  return *latest_tree_;
}

TreeCollection Driver::ParseString(const std::string &str) {
  Clear();
  yy::parser parser_instance(*this);
  parser_instance.set_debug_level(trace_parsing_);
  Tree::TreeVector trees = {ParseString(&parser_instance, str)};
  return TreeCollection(trees, this->TagTaxonMap());
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
