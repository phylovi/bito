// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SBN_TRAINING_HPP_
#define SRC_SBN_TRAINING_HPP_

#include "sbn_maps.hpp"

using IndexerRepresentationCounter =
    std::vector<std::pair<IndexerRepresentation, uint32_t>>;

namespace SBNTraining {

IndexerRepresentationCounter IndexerRepresentationCounterOf(
    const BitsetSizeMap& indexer, const Node::TopologyCounter& topology_counter);

// The "SBN-SA" estimator described in the "Maximum Lower Bound Estimates" section of
// the 2018 NeurIPS paper.
void SimpleAverage(std::vector<double>& sbn_parameters,
                   const IndexerRepresentationCounter& indexer_representation_counter);

}  // namespace SBNTraining

#ifdef DOCTEST_LIBRARY_INCLUDED



#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_SBN_TRAINING_HPP_
