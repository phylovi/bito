// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "sbn_maps.hpp"
#include "sbn_support.hpp"

class UnrootedSBNSupport : public SBNSupport {
 public:
  UnrootedSBNSupport() : SBNSupport({}){};
  explicit UnrootedSBNSupport(const Node::TopologyCounter &topologies,
                              StringVector taxon_names)
      : SBNSupport(std::move(taxon_names)) {
    std::tie(rootsplits_, indexer_, index_to_child_, parent_to_child_range_,
             gpcsp_count_) =
        SBNMaps::BuildIndexerBundle(UnrootedSBNMaps::RootsplitCounterOf(topologies),
                                    UnrootedSBNMaps::PCSPCounterOf(topologies));
  }

  UnrootedIndexerRepresentationCounter IndexerRepresentationCounterOf(
      const Node::TopologyCounter &topology_counter, const size_t out_of_sample_index) {
    return UnrootedSBNMaps::IndexerRepresentationCounterOf(indexer_, topology_counter,
                                                           out_of_sample_index);
  }

  UnrootedIndexerRepresentationCounter IndexerRepresentationCounterOf(
      const Node::TopologyCounter &topology_counter) {
    return IndexerRepresentationCounterOf(topology_counter, GPCSPCount());
  }

  UnrootedIndexerRepresentation IndexerRepresentationOf(
      const Node::NodePtr &topology, const size_t out_of_sample_index) const {
    return UnrootedSBNMaps::IndexerRepresentationOf(indexer_, topology,
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
