// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_UNROOTED_SBN_SUPPORT_HPP_
#define SRC_UNROOTED_SBN_SUPPORT_HPP_

#include "sbn_maps.hpp"
#include "sbn_support.hpp"

class UnrootedSBNSupport : public SBNSupport {
 public:
  UnrootedSBNSupport() : SBNSupport({}){};
  explicit UnrootedSBNSupport(const Node::TopologyCounter &topologies,
                              StringVector taxon_names)
      : SBNSupport(std::move(taxon_names)) {
    std::tie(rootsplits_, indexer_, index_to_child_, parent_to_range_, gpcsp_count_) =
        SBNMaps::BuildIndexerBundle(UnrootedSBNMaps::RootsplitSupportOf(topologies),
                                    UnrootedSBNMaps::PCSPSupportOf(topologies));
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

#endif  // SRC_UNROOTED_SBN_SUPPORT_HPP_
