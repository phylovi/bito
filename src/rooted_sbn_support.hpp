// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "sbn_support.hpp"

class RootedSBNSupport : public SBNSupport {
 public:
  RootedSBNSupport() : SBNSupport(nullptr) {};
  explicit RootedSBNSupport(SubsplitDAG* dag) : SBNSupport(dag) {}

  RootedIndexerRepresentationCounter IndexerRepresentationCounterOf(
      const Node::TopologyCounter &topology_counter, const size_t out_of_sample_index) {
    return RootedSBNMaps::IndexerRepresentationCounterOf(Indexer(), topology_counter,
                                                         out_of_sample_index);
  }

  RootedIndexerRepresentationCounter IndexerRepresentationCounterOf(
      const Node::TopologyCounter &topology_counter) {
    return IndexerRepresentationCounterOf(topology_counter, GPCSPCount());
  }

  RootedIndexerRepresentation IndexerRepresentationOf(
      const Node::NodePtr &topology, const size_t out_of_sample_index) const {
    return RootedSBNMaps::IndexerRepresentationOf(Indexer(), topology,
                                                  out_of_sample_index);
  }

  RootedIndexerRepresentation IndexerRepresentationOf(
      const Node::NodePtr &topology) const {
    return IndexerRepresentationOf(topology, GPCSPCount());
  }

  static BitsetSizeDict RootsplitCounterOf(const Node::TopologyCounter &topologies) {
    return RootedSBNMaps::RootsplitCounterOf(topologies);
  }

  static PCSPCounter PCSPCounterOf(const Node::TopologyCounter &topologies) {
    return RootedSBNMaps::PCSPCounterOf(topologies);
  }
};
