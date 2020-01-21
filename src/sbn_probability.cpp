// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "sbn_probability.hpp"
#include <iostream>
#include <numeric>
#include "ProgressBar.hpp"
#include "numerical_utils.hpp"
#include "sbn_maps.hpp"

// Increment all entries from an index vector by a log(value).
void IncrementByInLog(EigenVectorXdRef vec, const SizeVector& indices, double value) {
  for (const auto& idx : indices) {
    vec[idx] = NumericalUtils::LogAdd(vec[idx], value);
  }
}

// Increment all entries from an index vector vector by a log(value).
void IncrementByInLog(EigenVectorXdRef vec, const SizeVectorVector& index_vector_vector,
                      double value) {
  for (const auto& indices : index_vector_vector) {
    IncrementByInLog(vec, indices, value);
  }
}

// Increment all entries from an index vector by a value.
void IncrementBy(EigenVectorXdRef vec, const SizeVector& indices, double value) {
  for (const auto& idx : indices) {
    vec[idx] += value;
  }
}

// Increment all entries from an index vector vector by a value.
void IncrementBy(EigenVectorXdRef vec, const SizeVectorVector& index_vector_vector,
                 double value) {
  for (const auto& indices : index_vector_vector) {
    IncrementBy(vec, indices, value);
  }
}

// Increment all entries from an index vector of length k by a value taken from a vector
// of values of length k.
void IncrementBy(EigenVectorXdRef vec, const SizeVector& indices,
                 const EigenConstVectorXdRef values) {
  Assert(indices.size() == values.size(),
         "Indices and values don't have matching size.");
  for (size_t i = 0; i < values.size(); ++i) {
    vec[indices[i]] += values[i];
  }
}

// Repeat the previous increment operation across a vector of index vectors with a fixed
// vector of values.
void IncrementBy(EigenVectorXdRef vec, const SizeVectorVector& index_vector_vector,
                 const EigenConstVectorXdRef values) {
  Assert(index_vector_vector.size() == values.size(),
         "Indices and values don't have matching size.");
  for (size_t i = 0; i < values.size(); ++i) {
    IncrementBy(vec, index_vector_vector[i], values[i]);
  }
}

// Take the product of the entries of vec in indices times starting_value.
double ProductOf(const EigenConstVectorXdRef vec, const SizeVector& indices,
                 const double starting_value) {
  double result = starting_value;
  for (const auto& idx : indices) {
    result *= vec[idx];
  }
  return result;
}

double SumOf(const EigenConstVectorXdRef vec, const SizeVector& indices,
             const double starting_value) {
  double result = starting_value;
  for (const auto& idx : indices) {
    result += vec[idx];
  }
  return result;
}

// Probability-normalize a range of values in a vector.
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

void SBNProbability::ProbabilityNormalizeRangeInLog(EigenVectorXdRef vec,
                                                    std::pair<size_t, size_t> range) {
  auto [start_idx, end_idx] = range;
  auto segment = vec.segment(start_idx, end_idx - start_idx);
  NumericalUtils::ProbabilityNormalizeInLog(segment);
}

void SBNProbability::ProbabilityNormalizeParamsInLog(
    EigenVectorXdRef vec, size_t rootsplit_count,
    const BitsetSizePairMap& parent_to_range) {
  ProbabilityNormalizeRangeInLog(vec, {0, rootsplit_count});
  for (const auto& [_, range] : parent_to_range) {
    ProbabilityNormalizeRangeInLog(vec, range);
  }
}

// Set the provided counts vector to be the counts of the rootsplits and PCSSs provided
// in the input.
void SetCounts(EigenVectorXdRef counts,
               const IndexerRepresentationCounter& indexer_representation_counter,
               size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  counts.setZero();
  for (const auto& [indexer_representation, int_topology_count] :
       indexer_representation_counter) {
    const auto topology_count = static_cast<double>(int_topology_count);
    IncrementBy(counts, indexer_representation, topology_count);
  }
}

// Set the provided counts vector to be the counts of the rootsplits and PCSSs provided
// in the input.
void SetLogCounts(EigenVectorXdRef counts,
                  const IndexerRepresentationCounter& indexer_representation_counter,
                  size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  counts.fill(DOUBLE_NEG_INF);
  for (const auto& [indexer_representation, int_topology_count] :
       indexer_representation_counter) {
    const auto log_topology_count = log(static_cast<double>(int_topology_count));
    IncrementByInLog(counts, indexer_representation, log_topology_count);
  }
}

void SBNProbability::SimpleAverage(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  SetLogCounts(sbn_parameters, indexer_representation_counter, rootsplit_count,
               parent_to_range);
}

EigenVectorXd SBNProbability::ExpectationMaximization(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range, double alpha,
    size_t em_loop_count) {
  Assert(!indexer_representation_counter.empty(),
         "Empty indexer_representation_counter.");
  auto edge_count = indexer_representation_counter[0].first.size();
  // This vector holds the \bar{m} vectors (described in the 2018 NeurIPS paper).
  // They are packed into a single vector as sbn_parameters is.
  EigenVectorXd m_bar(sbn_parameters.size());
  // // The log of the sbn_parameters (only used for calculating the regularized score).
  // EigenVectorXd log_sbn_parameters(sbn_parameters.size());
  // The q weight of a rootsplit is the probability of each rooting given the current
  // SBN parameters.
  EigenVectorXd q_weights(edge_count);
  // This vector holds the \tilde{m} vectors (described in the 2018 NeurIPS paper),
  // which is the counts vector before normalization to get the SimpleAverage estimate.
  EigenVectorXd m_tilde(sbn_parameters.size());
  SetCounts(m_tilde, indexer_representation_counter, rootsplit_count, parent_to_range);
  // m_tilde is the counts, but marginalized over a uniform distribution on the rooting
  // edge. Thus we take the total counts and then divide by the edge count.
  m_tilde = m_tilde / static_cast<double>(edge_count);
  // The normalized version of m_tilde is the SA estimate, which is our starting point.
  sbn_parameters = m_tilde;
  ProbabilityNormalizeParams(sbn_parameters, rootsplit_count, parent_to_range);
  // The score is the Q^(n) defined on p.6 of the 2018 NeurIPS paper.
  EigenVectorXd score_history(em_loop_count);
  // Do the specified number of EM loops.
  ProgressBar progress_bar(em_loop_count, 70);
  for (size_t em_idx = 0; em_idx < em_loop_count; ++em_idx) {
    double score = 0.;
    m_bar.setZero();
    // Loop over topologies (as manifested by their indexer representations).
    for (const auto& [indexer_representation, int_topology_count] :
         indexer_representation_counter) {
      double unnormalized_partial_score = 0.;
      const auto topology_count = static_cast<double>(int_topology_count);
      // Calculate the q weights for this topology.
      q_weights.setZero();
      Assert(indexer_representation.size() == edge_count,
             "Indexer representation length is not constant.");
      // Loop over the various rooting positions of this topology.
      for (size_t rooting_position = 0; rooting_position < edge_count;
           ++rooting_position) {
        const SizeVector& rooted_representation =
            indexer_representation[rooting_position];
        // Calculate the SBN probability of this topology rooted at this position.
        double p_rooted_topology = ProductOf(sbn_parameters, rooted_representation, 1.);
        // The unnormalized q_weight is this probability.
        q_weights[rooting_position] = p_rooted_topology;
        // This is the score with an unnormalized q.
        unnormalized_partial_score += p_rooted_topology * log(p_rooted_topology);
      }  // End of looping over rooting positions.
      double q_sum = q_weights.sum();
      // Normalize q_weights.
      q_weights /= q_sum;
      // Now increment the new SBN parameters by the q-weighted counts.
      q_weights *= topology_count;
      IncrementBy(m_bar, indexer_representation, q_weights);
      // We increment the score with the per-topology-contribution after applying the
      // normalization factor for q, turning it into an actual probability.
      score += unnormalized_partial_score / q_sum;
    }  // End of looping over topologies.
    sbn_parameters = m_bar + alpha * m_tilde;
    ProbabilityNormalizeParams(sbn_parameters, rootsplit_count, parent_to_range);
    // if (alpha > 0.) {
    //   log_sbn_parameters = sbn_parameters.log();
    //   score += alpha * m_tilde.dot(log_sbn_parameters);
    // }
    score_history[em_idx] = score;
    ++progress_bar;
    progress_bar.display();
  }  // End of EM loop.
  progress_bar.done();
  sbn_parameters = sbn_parameters.array().log();
  return score_history;
}

double SBNProbability::ProbabilityOf(
    const EigenConstVectorXdRef sbn_parameters,
    const IndexerRepresentation& indexer_representation) {
  double log_total_probability = DOUBLE_NEG_INF;
  for (const auto& rooted_representation : indexer_representation) {
    log_total_probability = NumericalUtils::LogAdd(
        log_total_probability, SumOf(sbn_parameters, rooted_representation, 0.));
  };
  return exp(log_total_probability);
}

EigenVectorXd SBNProbability::ProbabilityOf(
    const EigenConstVectorXdRef sbn_parameters,
    const std::vector<IndexerRepresentation>& indexer_representations) {
  const size_t topology_count = indexer_representations.size();
  EigenVectorXd results(topology_count);
  for (size_t topology_idx = 0; topology_idx < topology_count; ++topology_idx) {
    results[topology_idx] =
        ProbabilityOf(sbn_parameters, indexer_representations[topology_idx]);
  }
  return results;
}
