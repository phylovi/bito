// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A simple version of the SBN support, which will be used to build SBNSupport objects.
// Here we have no indexer functionality, and no taxon names.
//
// Factoring things out in this way allows us to have simple tests for parsing.

#ifndef SRC_SIMPLE_SBN_SUPPORT_HPP_
#define SRC_SIMPLE_SBN_SUPPORT_HPP_

#include "sbn_maps.hpp"

class SimpleSBNSupport {
 public:
  SimpleSBNSupport(BitsetSet rootsplit_support, PCSPSupportDict pcsp_support)
      : rootsplit_support_(std::move(rootsplit_support)),
        pcsp_support_(std::move(pcsp_support)){};
  SimpleSBNSupport(BitsetSizeDict rootsplit_counter, PCSPCounter pcsp_counter);

  IndexerBundle ToIndexerBundle();

  // static OfRootedTopologies
  //     RootedSBNMaps::RootsplitCounterOf(topologies),
  //     RootedSBNMaps::PCSPCounterOf(topologies)
  // TODO NEXT: convert BuildIndexerBundle to take this ensemble of things, then convert
  // to take these static functions.
  //
  // static OfUnrootedTopologies
  // static OfPrettyStringDoubleMap
  // static OfPrettySBNParameterCSV

 protected:
  BitsetSet rootsplit_support_;
  PCSPSupportDict pcsp_support_;
};

#endif  // SRC_SIMPLE_SBN_SUPPORT_HPP_
