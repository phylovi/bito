// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "sbn_training.hpp"
#include "sbn_maps.hpp"

void IncrementBy(EigenVectorXdRef vec, const SizeVector& indices, double value) {
  for (const auto& idx : indices) {
    vec[idx] += value;
  }
}

void IncrementBy(EigenVectorXdRef vec, const SizeVectorVector& index_vector_vector,
                 double value) {
  for (const auto& indices : index_vector_vector) {
    IncrementBy(vec, indices, value);
  }
}

IndexerRepresentationCounter SBNTraining::IndexerRepresentationCounterOf(
    const BitsetSizeMap& indexer, const Node::TopologyCounter& topology_counter) {
  IndexerRepresentationCounter counter;
  counter.reserve(topology_counter.size());
  for (const auto& [topology, count] : topology_counter) {
    counter.push_back({SBNMaps::IndexerRepresentationOf(indexer, topology), count});
  }
  return counter;
}

void SBNTraining::SimpleAverage(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter) {
  sbn_parameters.setZero();

  for (const auto& [indexer_representation, int_count] :
       indexer_representation_counter) {
    const auto& [rootsplits, pcsss] = indexer_representation;
    const auto count = static_cast<double>(int_count);

    IncrementBy(sbn_parameters, rootsplits, count);
    IncrementBy(sbn_parameters, pcsss, count);
  }
}
