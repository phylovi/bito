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

