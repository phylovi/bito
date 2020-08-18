// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_SBN_SUPPORT_HPP_
#define SRC_ROOTED_SBN_SUPPORT_HPP_

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

  static BitsetSizeDict RootsplitCounterOf(const Node::TopologyCounter &topologies) {
    return RootedSBNMaps::RootsplitCounterOf(topologies);
  }
  static PCSPDict PCSPCounterOf(const Node::TopologyCounter &topologies) {
    return RootedSBNMaps::PCSPCounterOf(topologies);
  }
};

#endif  // SRC_ROOTED_SBN_SUPPORT_HPP_
