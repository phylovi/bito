// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "sbn_maps.hpp"
#include "sbn_support.hpp"

class UnrootedSBNSupport : public SBNSupport {
 public:
  UnrootedSBNSupport() : SBNSupport(nullptr) {};
  explicit UnrootedSBNSupport(SubsplitDAG* dag) : SBNSupport(dag) {}

  UnrootedIndexerRepresentationCounter IndexerRepresentationCounterOf(
      const Node::TopologyCounter &topology_counter, const size_t out_of_sample_index) {
    return UnrootedSBNMaps::IndexerRepresentationCounterOf(Indexer(), topology_counter,
                                                           out_of_sample_index);
  }

  UnrootedIndexerRepresentationCounter IndexerRepresentationCounterOf(
      const Node::TopologyCounter &topology_counter) {
    return IndexerRepresentationCounterOf(topology_counter, GPCSPCount());
  }

  UnrootedIndexerRepresentation IndexerRepresentationOf(
      const Node::NodePtr &topology, const size_t out_of_sample_index) const {
    return UnrootedSBNMaps::IndexerRepresentationOf(Indexer(), topology,
                                                    out_of_sample_index);
  }

  UnrootedIndexerRepresentation IndexerRepresentationOf(
      const Node::NodePtr &topology) const {
    return IndexerRepresentationOf(topology, GPCSPCount());
  }

  static BitsetSizeDict RootsplitCounterOf(const Node::TopologyCounter &topologies) {
    return UnrootedSBNMaps::RootsplitCounterOf(topologies);
  }

  static PCSPCounter PCSPCounterOf(const Node::TopologyCounter &topologies) {
    return UnrootedSBNMaps::PCSPCounterOf(topologies);
  }
};
