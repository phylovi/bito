// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Perform training of an SBN based on a sample of trees.
//
// We assume that readers are familiar with how the sbn_parameters_ vector is laid out:
// first probabilities of rootsplits, then conditional probabilities of PCSSs.

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
  for (size_t i = 0; i < values.size(); ++i) {
    vec[indices[i]] += values[i];
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
  double result = 1.;
  for (const auto& idx : indices) {
    result *= vec[idx];
  }
  return result;
}

void ProbabilityNormalizeRange(EigenVectorXdRef vec, std::pair<size_t, size_t> range) {
  auto [start_idx, end_idx] = range;
  auto segment = vec.segment(start_idx, end_idx - start_idx);
  segment /= segment.sum();
}

// We assume that vec is laid out like sbn_parameters (see top).
void ProbabilityNormalizeParams(EigenVectorXdRef vec, size_t rootsplit_count,
                                const BitsetSizePairMap& parent_to_range) {
  ProbabilityNormalizeRange(vec, {0, rootsplit_count});
  for (const auto& [_, range] : parent_to_range) {
    ProbabilityNormalizeRange(vec, range);
  }
}

IndexerRepresentationCounter SBNTraining::IndexerRepresentationCounterOf(
    const BitsetSizeMap& indexer, const Node::TopologyCounter& topology_counter) {
  IndexerRepresentationCounter counter;
  counter.reserve(topology_counter.size());
  for (const auto& [topology, topology_count] : topology_counter) {
    counter.push_back(
        {SBNMaps::IndexerRepresentationOf(indexer, topology), topology_count});
  }
  return counter;
}

void SBNTraining::SimpleAverage(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter) {
  sbn_parameters.setZero();

  for (const auto& [indexer_representation, int_topology_count] :
       indexer_representation_counter) {
    const auto& [rootsplits, pcsss] = indexer_representation;
    const auto topology_count = static_cast<double>(int_topology_count);
    IncrementBy(sbn_parameters, rootsplits, topology_count);
    IncrementBy(sbn_parameters, pcsss, topology_count);
  }
  // We leave the counts in an un-normalized state, which is a suitable input for
  // sampling.
  // TODO but it's not suitable for probability normalization.
}

void SBNTraining::ExpectationMaximization(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range, double alpha,
    double tolerance) {
  // TODO use alpha and tolerance
  // This vector holds the \bar{m} vectors (described in the 2018 NeurIPS paper).
  // They are packed into a single vector as sbn_parameters is.
  EigenVectorXd m_bar(sbn_parameters.size());
  // The q weight of a rootsplit is the probability of each rooting given the current
  // SBN parameters.
  Assert(!indexer_representation_counter.empty(),
         "Empty indexer_representation_counter.");
  auto edge_count = indexer_representation_counter[0].first.first.size();
  EigenVectorXd q_weights(edge_count);

  SimpleAverage(sbn_parameters, indexer_representation_counter);
  ProbabilityNormalizeParams(sbn_parameters, rootsplit_count, parent_to_range);

  // TODO actually loop
  // Start EM loop.
  m_bar.setZero();
  // Loop over topologies (as manifested by their indexer representations).
  for (const auto& [indexer_representation, int_topology_count] :
       indexer_representation_counter) {
    const auto& [rootsplits, pcsss] = indexer_representation;
    const auto topology_count = static_cast<double>(int_topology_count);
    // Calculate the q weights for this topology.
    q_weights.setZero();
    Assert(rootsplits.size() == edge_count,
           "Rootsplit length not equal to edge_count.");
    Assert(pcsss.size() == edge_count, "PCSSs length not equal to edge_count.");
    // Loop over the various rooting positions of this topology.
    for (size_t rooting_position = 0; rooting_position < edge_count;
         ++rooting_position) {
      const size_t& rootsplit = rootsplits[rooting_position];
      const SizeVector& pcss = pcsss[rooting_position];
      // Calculate the SBN probability of this topology rooted at this position.
      q_weights[rooting_position] =
          sbn_parameters[rootsplit] * ProductOf(sbn_parameters, pcss);
    }
    q_weights /= q_weights.sum();
    // Now increment the new SBN parameters by the q-weighted counts.
    q_weights *= topology_count;
    IncrementBy(m_bar, rootsplits, q_weights);
    IncrementBy(m_bar, pcsss, q_weights);
  }
  sbn_parameters = m_bar;
  ProbabilityNormalizeParams(sbn_parameters, rootsplit_count, parent_to_range);
}
