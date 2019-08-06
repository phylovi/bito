// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "alignment.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>

size_t Alignment::Length() const {
  Assert(SequenceCount() > 0,
         "Must have sequences in an alignment to ask for a Length.");
  return data_.begin()->second.size();
}

bool Alignment::IsValid() const {
  if (data_.size() == 0) {
    return false;
  }
  size_t length = data_.begin()->second.size();
  for (const auto &iter : data_) {
    if (length != iter.second.size()) {
      return false;
    }
  }
  return true;
}

std::string Alignment::at(const std::string &taxon) const {
  auto search = data_.find(taxon);
  if (search != data_.end()) {
    return search->second;
  } else {
    Failwith("Taxon '" + taxon + "' not found in alignment.");
  }
}

// An edited version of
// https://stackoverflow.com/questions/35251635/fasta-reader-written-in-c
// which seems like it was originally taken from
// http://rosettacode.org/wiki/FASTA_format#C.2B.2B
void Alignment::ReadFasta(std::string fname) {
  StringStringMap &data = this->data_;
  data.clear();
  auto insert = [&data](std::string taxon, std::string sequence) {
    if (!taxon.empty()) {
      SafeInsert(data, taxon, sequence);
    }
  };
  std::ifstream input(fname);
  if (!input.good()) {
    Failwith("Could not open '" + fname + "'");
  }
  std::string line, taxon, sequence;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    if (line[0] == '>') {
      insert(taxon, sequence);
      taxon = line.substr(1);
      sequence.clear();
    } else {
      sequence += line;
    }
  }
  // Insert the last taxon, sequence pair.
  insert(taxon, sequence);
  if (!IsValid()) {
    Failwith("Sequences of the alignment are not all the same length.");
  }
}
