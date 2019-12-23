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

void IncrementBy(EigenVectorXdRef vec, const SizeVector& indices,
                 const EigenConstVectorXdRef values) {
  Assert(indices.size() == values.size(),
         "Indices and values don't have matching size.");
  for (const auto& idx : indices) {
    vec[idx] += values[idx];
  }
}

void IncrementBy(EigenVectorXdRef vec, const SizeVectorVector& index_vector_vector,
                 const EigenConstVectorXdRef values) {
  Assert(index_vector_vector.size() == values.size(),
         "Indices and values don't have matching size.");
  for (size_t i = 0; i < values.size(); ++i) {
    IncrementBy(vec, index_vector_vector[i], values[i]);
  }
}

double ProductOf(EigenVectorXdRef vec, const SizeVector& indices) {
  double result = 0.;
  for (const auto& idx : indices) {
    result *= vec[idx];
  }
  return result;
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

void SBNTraining::ExpectationMaximization(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, double tolerance) {
  EigenVectorXd new_sbn_parameters(sbn_parameters.size());
  // The q weight of a rootsplit is the probability of each rooting given the current
  // SBN parameters.
  Assert(!indexer_representation_counter.empty(),
         "Empty indexer_representation_counter.");
  auto edge_count = indexer_representation_counter[0].first.first.size();
  EigenVectorXd q_weights(edge_count);

  SimpleAverage(sbn_parameters, indexer_representation_counter);
  // TODO normalize the rootsplit and SBN probabilities.

  // Start EM loop.
  new_sbn_parameters.setZero();
  // Loop over topologies (as manifested by their indexer representations).
  for (const auto& [indexer_representation, int_count] :
       indexer_representation_counter) {
    const auto& [rootsplits, pcsss] = indexer_representation;
    const auto count = static_cast<double>(int_count);
    // Calculate the q weights for this topology.
    q_weights.setZero();
    Assert(rootsplits.size() == pcsss.size(),
           "Rootsplit length not the same as pcss length.");
    // Loop over the various rooting positions of this topology.
    for (size_t rooting_position = 0; rooting_position < rootsplits.size();
         ++rooting_position) {
      const auto& rootsplit = rootsplits[rooting_position];
      const auto& pcss = pcsss[rooting_position];
      // Calculate the SBN probability of this topology rooted at this position.
      q_weights[rooting_position] =
          sbn_parameters[rootsplit] * ProductOf(sbn_parameters, pcss);
    }
    q_weights /= q_weights.sum();
    // Now increment the new SBN parameters by the q-weighted counts.
    q_weights *= count;
    IncrementBy(new_sbn_parameters, rootsplits, q_weights);
    IncrementBy(new_sbn_parameters, pcsss, q_weights);
  }
  sbn_parameters = new_sbn_parameters;
}
