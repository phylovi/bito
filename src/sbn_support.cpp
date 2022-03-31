// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "sbn_support.hpp"

StringVector SBNSupport::PrettyIndexer() const {
  auto indexer = Indexer();
  StringVector pretty_representation(indexer.size());
  for (const auto& [key, idx] : indexer) {
    pretty_representation[idx] = key.PCSPToString();
  }
  return pretty_representation;
}

void SBNSupport::PrettyPrintIndexer() const {
  auto pretty_representation = PrettyIndexer();
  for (size_t i = 0; i < pretty_representation.size(); i++) {
    std::cout << i << "\t" << pretty_representation[i] << std::endl;
  }
}

// Return indexer_ and parent_to_child_range_ converted into string-keyed maps.
std::tuple<StringSizeMap, StringSizePairMap> SBNSupport::GetIndexers() const {
  auto indexer = Indexer();
  auto str_indexer = StringifyMap(indexer);
  auto str_parent_to_range = StringifyMap(ParentToRange());
  std::string rootsplit("DAG Root Node");
  SafeInsert(str_parent_to_range, rootsplit, {0, RootsplitCount()});
  return {str_indexer, str_parent_to_range};
}

// Get the indexer, but reversed and with bitsets appropriately converted to
// strings.
StringVector SBNSupport::StringReversedIndexer() const {
  auto indexer = Indexer();
  StringVector reversed_indexer(indexer.size());
  for (const auto& [key, idx] : indexer) {
    reversed_indexer[idx] = key.PCSPToString();
  }
  return reversed_indexer;
}

void SBNSupport::ProbabilityNormalizeSBNParametersInLog(
    EigenVectorXdRef sbn_parameters) const {
  SBNProbability::ProbabilityNormalizeParamsInLog(sbn_parameters, RootsplitCount(),
                                                  ParentToRange());
}
