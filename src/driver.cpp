// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

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
#include "taxon_name_munging.hpp"
#include "zlib_stream.hpp"

Driver::Driver()
    : next_id_(0),
      sort_taxa_(false),
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
TreeCollection Driver::ParseNewick(std::istream &in) {
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
      // If taxon map has not be initialized, we parse the first tree temporarily to get
      // the taxon names.  We will re-parse this tree after assigning sorted taxon IDs.
      if (!taxa_complete_ && sort_taxa_) {
        ParseString(&parser_instance, line);
        SortTaxa();
      }
      trees.push_back(ParseString(&parser_instance, line));
    }
  }
  return TreeCollection(std::move(trees), this->TagTaxonMap());
}

TreeCollection Driver::ParseAndDequoteNewick(std::istream &in) {
  TreeCollection perhaps_quoted_trees = ParseNewick(in);
  return TreeCollection(
      std::move(perhaps_quoted_trees.trees_),
      TaxonNameMunging::DequoteTagStringMap(perhaps_quoted_trees.TagTaxonMap()));
}

TreeCollection Driver::ParseNewickFile(const std::string &fname) {
  Clear();
  std::ifstream in(fname.c_str());
  if (!in) {
    Failwith("Cannot open the File : " + fname);
  }
  return ParseAndDequoteNewick(in);
}

TreeCollection Driver::ParseNewickFileGZ(const std::string &fname) {
  Clear();
  std::ifstream in_compressed(fname.c_str());
  if (!in_compressed) {
    Failwith("Cannot open the File : " + fname);
  }
  zlib::ZStringBuf zbuf(in_compressed, 1024, 2048);
  std::istream in(&zbuf);
  return ParseAndDequoteNewick(in);
}

void GetLineAndConvertToLowerCase(std::istream &in, std::string &line) {
  std::getline(in, line);
  std::transform(line.begin(), line.end(), line.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}

TreeCollection Driver::ParseNexusFile(const std::string &fname) {
  std::ifstream in(fname.c_str());
  if (!in) {
    throw std::runtime_error("Cannot open file.");
  }
  return ParseNexus(in);
}

TreeCollection Driver::ParseNexusFileGZ(const std::string &fname) {
  std::ifstream in_compressed(fname.c_str());
  if (!in_compressed) {
    throw std::runtime_error("Cannot open file.");
  }
  zlib::ZStringBuf zbuf(in_compressed, 1024, 2048);
  std::istream in(&zbuf);
  return ParseNexus(in);
}

TreeCollection Driver::ParseNexus(std::istream &in) {
  Clear();
  try {
    std::string line;
    std::getline(in, line);
    if (line != "#NEXUS") {
      throw std::runtime_error("Putative Nexus file doesn't begin with #NEXUS.");
    }
    do {
      if (in.eof()) {
        throw std::runtime_error("Finished reading and couldn't find 'begin trees;'");
      }
      GetLineAndConvertToLowerCase(in, line);
    } while (line != "begin trees;");
    GetLineAndConvertToLowerCase(in, line);
    std::regex translate_start("^\\s*translate");
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
    // Iterate through the translate table, assigning tags according to the order of
    // taxa in the block. So, the first taxon name gets leaf number 0, etc.
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
                          TaxonNameMunging::DequoteTagStringMap(long_name_taxon_map));
  } catch (const std::exception &exception) {
    Failwith(std::string("Problem parsing Nexus file:\n") + exception.what());
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
  return TreeCollection(trees,
                        TaxonNameMunging::DequoteTagStringMap(this->TagTaxonMap()));
}

void Driver::SetTaxa(const std::map<std::string, uint32_t> taxa) {
  taxa_ = taxa;
  taxa_complete_ = true;
}

void Driver::InitializeTaxa() {
  std::map<std::string, uint32_t> taxa;
  uint32_t new_id = 0;
  for (const auto &[name, old_id] : taxa_) {
    taxa[name] = new_id;
    new_id++;
  }
  taxa_ = taxa;
  taxa_complete_ = true;
}

TagStringMap Driver::TagTaxonMap() {
  TagStringMap m;
  for (const auto &iter : taxa_) {
    // These are leaves, so the number of leaves below is 1.
    m[PackInts(iter.second, 1)] = iter.first;
  }
  return m;
}

void Driver::SetTaxa(const std::map<std::string, uint32_t> taxa) {
  taxa_ = taxa;
  taxa_complete_ = true;
}

void Driver::SortTaxa() {
  std::map<std::string, uint32_t> taxa;
  uint32_t new_id = 0;
  for (const auto &[name, old_id] : taxa_) {
    std::ignore = old_id;
    taxa[name] = new_id;
    new_id++;
  }
  SetTaxa(taxa);
}

// Note that a number of Driver methods are implemented in scanner.ll.
