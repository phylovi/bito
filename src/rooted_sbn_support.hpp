// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "sbn_support.hpp"

class RootedSBNSupport : public SBNSupport {
 public:
  RootedSBNSupport() : SBNSupport({}){};
  explicit RootedSBNSupport(const Node::TopologyCounter &topologies,
                            StringVector taxon_names)
      : SBNSupport(std::move(taxon_names)) {
    std::tie(rootsplits_, indexer_, index_to_child_, parent_to_range_, gpcsp_count_) =
        SBNMaps::BuildIndexerBundle(RootedSBNMaps::RootsplitCounterOf(topologies),
                                    RootedSBNMaps::PCSPCounterOf(topologies));
  }

  RootedIndexerRepresentationCounter IndexerRepresentationCounterOf(
      const Node::TopologyCounter &topology_counter, const size_t out_of_sample_index) {
    return RootedSBNMaps::IndexerRepresentationCounterOf(indexer_, topology_counter,
                                                         out_of_sample_index);
  }

  RootedIndexerRepresentationCounter IndexerRepresentationCounterOf(
      const Node::TopologyCounter &topology_counter) {
    return IndexerRepresentationCounterOf(topology_counter, GPCSPCount());
  }

  RootedIndexerRepresentation IndexerRepresentationOf(
      const Node::NodePtr &topology, const size_t out_of_sample_index) const {
    return RootedSBNMaps::IndexerRepresentationOf(indexer_, topology,
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

