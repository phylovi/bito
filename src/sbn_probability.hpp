// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SBN_PROBABILITY_HPP_
#define SRC_SBN_PROBABILITY_HPP_

#include "eigen_sugar.hpp"
#include "sbn_maps.hpp"

using IndexerRepresentationCounter =
    std::vector<std::pair<IndexerRepresentation, uint32_t>>;

namespace SBNProbability {

// Turn a TopologyCounter into an IndexerRepresentationCounter.
IndexerRepresentationCounter IndexerRepresentationCounterOf(
    const BitsetSizeMap& indexer, const Node::TopologyCounter& topology_counter);

// The "SBN-SA" estimator described in the "Maximum Lower Bound Estimates" section of
// the 2018 NeurIPS paper.
void SimpleAverage(EigenVectorXdRef sbn_parameters,
                   const IndexerRepresentationCounter& indexer_representation_counter,
                   size_t rootsplit_count, const BitsetSizePairMap& parent_to_range);

// The "SBN-EM" estimator described in the "Expectation Maximization" section of
// the 2018 NeurIPS paper.
void ExpectationMaximization(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range, double alpha,
    size_t em_loop_count);

// Calculate the probability of an indexer_representation of a topology.
double ProbabilityOf(const EigenConstVectorXdRef,
                     const IndexerRepresentation& indexer_representation);

// Calculate the probabilities of a collection of indexer_representations.
EigenVectorXd ProbabilityOf(
    const EigenConstVectorXdRef sbn_parameters,
    const std::vector<IndexerRepresentation>& indexer_representations);

}  // namespace SBNProbability

#ifdef DOCTEST_LIBRARY_INCLUDED



#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_SBN_PROBABILITY_HPP_
