// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "sbn_training.hpp"
#include "sbn_maps.hpp"

IndexerRepresentationCounter SBNTraining::IndexerRepresentationCounterOf(
    const BitsetSizeMap& indexer, const Node::TopologyCounter& topology_counter) {
  IndexerRepresentationCounter counter;
  counter.reserve(topology_counter.size());
  for (const auto& [topology, count] : topology_counter) {
    counter.push_back({SBNMaps::IndexerRepresentationOf(indexer, topology), count});
  }
  return counter;
}

void SBNTraining::SimpleAverage(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter) {
  sbn_parameters.setZero();

  for (const auto& [indexer_representation, int_count] :
       indexer_representation_counter) {
    const auto& [rootsplits, pcsss] = indexer_representation;
    const auto count = static_cast<double>(int_count);

    for (const auto& rootsplit_index : rootsplits) {
      sbn_parameters[rootsplit_index] += count;
    }
    for (const auto& pcss : pcsss) {
      for (const auto& pcss_index : pcss) {
        sbn_parameters[pcss_index] += count;
      }
    }
  }
}
