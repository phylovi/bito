// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

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
  if (data_.empty()) {
    return false;
  }
  auto length = Length();
  auto is_same_length = [length](const auto &datum) {
    return (length == datum.second.size());
  };
  return std::all_of(data_.cbegin(), data_.cend(), is_same_length);
}

const std::string &Alignment::at(const std::string &taxon) const {
  auto search = data_.find(taxon);
  if (search != data_.end()) {
    return search->second;
  }  // else
  Failwith("Taxon '" + taxon + "' not found in alignment.");
}

// An edited version of
// https://stackoverflow.com/questions/35251635/fasta-reader-written-in-c
// which seems like it was originally taken from
// http://rosettacode.org/wiki/FASTA_format#C.2B.2B
Alignment Alignment::ReadFasta(const std::string &fname) {
  StringStringMap data;
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
    if (line.empty()) {
      continue;
    }
    // else:
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
  Alignment alignment(data);
  if (!alignment.IsValid()) {
    Failwith("Sequences of the alignment are not all the same length.");
  }
  return alignment;
}

Alignment Alignment::ExtractSingleColumnAlignment(size_t which_column) const {
  Assert(which_column < Length(),
         "Alignment::ExtractSingleColumnAlignment: Given column is longer than "
         "sequence length.");
  StringStringMap out_map;
  for (const auto &[taxon_name, sequence] : data_) {
    SafeInsert(out_map, taxon_name, sequence.substr(which_column, 1));
  }
  return Alignment(out_map);
}
