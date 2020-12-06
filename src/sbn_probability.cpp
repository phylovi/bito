// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "sbn_probability.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

#include "ProgressBar.hpp"
#include "numerical_utils.hpp"

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

// Increment all entries from an index vector of length k by a value taken from a vector
// of values of length k.
void IncrementByInLog(EigenVectorXdRef vec, const SizeVector& indices,
                      const EigenConstVectorXdRef values) {
  Assert(indices.size() == values.size(),
         "Indices and values don't have matching size.");
  for (size_t i = 0; i < values.size(); ++i) {
    vec[indices[i]] = NumericalUtils::LogAdd(vec[indices[i]], values[i]);
  }
}

// Repeat the previous increment operation across a vector of index vectors with a fixed
// vector of values.
void IncrementByInLog(EigenVectorXdRef vec, const SizeVectorVector& index_vector_vector,
                      const EigenConstVectorXdRef values) {
  Assert(index_vector_vector.size() == values.size(),
         "Indices and values don't have matching size.");
  for (size_t i = 0; i < values.size(); ++i) {
    IncrementByInLog(vec, index_vector_vector[i], values[i]);
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

// Take the sum of the entries of vec in indices plus starting_value.
double SBNProbability::SumOf(const EigenConstVectorXdRef vec, const SizeVector& indices,
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

// Normalize such that values become logs of probabilities.
void SBNProbability::ProbabilityNormalizeRangeInLog(EigenVectorXdRef vec,
                                                    std::pair<size_t, size_t> range) {
  auto [start_idx, end_idx] = range;
  auto segment = vec.segment(start_idx, end_idx - start_idx);
  NumericalUtils::ProbabilityNormalizeInLog(segment);
}

// We assume that vec is laid out like sbn_parameters (see top).
void SBNProbability::ProbabilityNormalizeParamsInLog(
    EigenVectorXdRef vec, size_t rootsplit_count,
    const BitsetSizePairMap& parent_to_range) {
  ProbabilityNormalizeRangeInLog(vec, {0, rootsplit_count});
  for (const auto& [_, range] : parent_to_range) {
    ProbabilityNormalizeRangeInLog(vec, range);
  }
}

// Set the provided counts vector to be the counts of the rootsplits and PCSPs provided
// in the input.
void SetCounts(
    EigenVectorXdRef counts,
    const UnrootedIndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  counts.setZero();
  for (const auto& [indexer_representation, int_topology_count] :
       indexer_representation_counter) {
    const auto topology_count = static_cast<double>(int_topology_count);
    IncrementBy(counts, indexer_representation, topology_count);
  }
}

// Set the provided counts vector to be the log of the counts of the rootsplits and
// PCSPs provided in the input.
// Note code duplication with the override below for
// UnrootedIndexerRepresentationCounter.
// We could refactor with templates, but then this whole file would have to go in the
// header.
void SetLogCounts(
    EigenVectorXdRef counts,
    const RootedIndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  counts.fill(DOUBLE_NEG_INF);
  for (const auto& [indexer_representation, int_topology_count] :
       indexer_representation_counter) {
    const auto log_topology_count = log(static_cast<double>(int_topology_count));
    IncrementByInLog(counts, indexer_representation, log_topology_count);
  }
}

// Note code duplication with the override below for
// UnrootedIndexerRepresentationCounter.
void SBNProbability::SimpleAverage(
    EigenVectorXdRef sbn_parameters,
    const RootedIndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  SetLogCounts(sbn_parameters, indexer_representation_counter, rootsplit_count,
               parent_to_range);
}

// Set the provided counts vector to be the log of the counts of the rootsplits and
// PCSPs provided in the input.
// Note code duplication with the override above for RootedIndexerRepresentationCounter.
void SetLogCounts(
    EigenVectorXdRef counts,
    const UnrootedIndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  counts.fill(DOUBLE_NEG_INF);
  for (const auto& [indexer_representation, int_topology_count] :
       indexer_representation_counter) {
    const auto log_topology_count = log(static_cast<double>(int_topology_count));
    IncrementByInLog(counts, indexer_representation, log_topology_count);
  }
}

// Note code duplication with the override above for RootedIndexerRepresentationCounter.
void SBNProbability::SimpleAverage(
    EigenVectorXdRef sbn_parameters,
    const UnrootedIndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  SetLogCounts(sbn_parameters, indexer_representation_counter, rootsplit_count,
               parent_to_range);
}

// All references to equations, etc, are to the 2018 NeurIPS paper.
// However, if you are doing a detailed read see doc/tex, because our definition of
// score differs from that in the NeurIPS paper, and also for details of how the prior
// calculation works.
EigenVectorXd SBNProbability::ExpectationMaximization(
    EigenVectorXdRef sbn_parameters,
    const UnrootedIndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range, double alpha,
    size_t max_iter, double score_epsilon) {
  Assert(!indexer_representation_counter.empty(),
         "Empty indexer_representation_counter.");
  auto edge_count = indexer_representation_counter[0].first.size();
  // The \bar{m} vectors (Algorithm 1) in log space.
  // They are packed into a single vector as sbn_parameters is.
  EigenVectorXd log_m_bar(sbn_parameters.size());
  // The q weight of a rootsplit is the probability of each rooting given the current
  // SBN parameters.
  EigenVectorXd log_q_weights(edge_count);
  // The \tilde{m} vectors (p.6): the counts vector before normalization to get the
  // SimpleAverage estimate. If alpha is nonzero log_m_tilde gets scaled by it below.
  EigenVectorXd log_m_tilde(sbn_parameters.size());
  SetLogCounts(log_m_tilde, indexer_representation_counter, rootsplit_count,
               parent_to_range);
  // m_tilde is the counts, but marginalized over a uniform distribution on the rooting
  // edge. Thus we take the total counts and then divide by the edge count.
  log_m_tilde = log_m_tilde.array() - log(static_cast<double>(edge_count));
  // The normalized version of m_tilde is the SA estimate, which is our starting point.
  sbn_parameters = log_m_tilde;
  // We need to ensure sbn_parameters is normalized as we are computing log P(S_1, T^u)
  // repeatedly.
  ProbabilityNormalizeParamsInLog(sbn_parameters, rootsplit_count, parent_to_range);
  // We need an exponentiated version of log_m_tilde for the score calculation if alpha
  // is nonzero.
  EigenVectorXd m_tilde_for_positive_alpha;
  if (alpha > 0.) {
    // For the regularized case, we always need log(alpha) + log_m_tilde so we store
    // this in log_m_tilde.
    log_m_tilde = log_m_tilde.array() + log(alpha);
    // We also need exp(log_m_tilde) = \alpha * tilde{m}_{s|t} for the regularized EM
    // algorithm.
    m_tilde_for_positive_alpha = log_m_tilde.array().exp();
  }
  // Our score is the marginal log likelihood of the training collection of trees (see
  // doc/tex).
  EigenVectorXd score_history = EigenVectorXd::Zero(max_iter);
  // Do the specified number of EM loops.
  ProgressBar progress_bar(max_iter);
  for (size_t em_idx = 0; em_idx < max_iter; ++em_idx) {
    log_m_bar.setConstant(DOUBLE_NEG_INF);
    // Loop over topologies (as manifested by their indexer representations).
    for (const auto& [indexer_representation, int_topology_count] :
         indexer_representation_counter) {
      // The number of times this topology was seen in the counter.
      const auto topology_count = static_cast<double>(int_topology_count);
      // Calculate the q weights for this topology.
      log_q_weights.setConstant(DOUBLE_NEG_INF);
      Assert(indexer_representation.size() == edge_count,
             "Indexer representation length is not constant.");
      // Loop over the various rooting positions of this topology, using log_q_weights
      // to store the probability of the tree in the various rootings (we will normalize
      // it later).
      for (size_t rooting_position = 0; rooting_position < edge_count;
           ++rooting_position) {
        const RootedIndexerRepresentation& rooted_representation =
            indexer_representation[rooting_position];
        // Calculate the SBN probability of this topology rooted at this position.
        double log_p_rooted_topology = SumOf(sbn_parameters, rooted_representation, 0.);
        // SHJ: Sometimes overflow is reported, sometimes it's underflow...
        if (fetestexcept(FE_OVER_AND_UNDER_FLOW_EXCEPT)) {
          log_q_weights[rooting_position] = DOUBLE_MINIMUM;
          feclearexcept(FE_OVER_AND_UNDER_FLOW_EXCEPT);
        } else {
          log_q_weights[rooting_position] = log_p_rooted_topology;
        }
      }  // End of looping over rooting positions.
      double log_p_unrooted_topology = NumericalUtils::LogSum(log_q_weights);
      score_history[em_idx] += topology_count * log_p_unrooted_topology;
      // Normalize q_weights to achieve the E-step of Algorithm 1.
      // For the increment step (M-step of Algorithm 1) we want a full topology
      // count rather than just the unique count. So we multiply the q_weights by the
      // topology count (in log space, it becomes summation rather than multiplication).
      log_q_weights =
          log_q_weights.array() + (-log_p_unrooted_topology + log(topology_count));
      // Increment the SBN-parameters-to-be by the q-weighted counts.
      IncrementByInLog(log_m_bar, indexer_representation, log_q_weights);
    }  // End of looping over topologies.
    // Store the proper value in sbn_parameters.
    sbn_parameters = (alpha > 0.)
                         ? NumericalUtils::LogAddVectors(log_m_bar, log_m_tilde)
                         : log_m_bar;
    // We normalize sbn_parameters right away to ensure that it is always normalized.
    ProbabilityNormalizeParamsInLog(sbn_parameters, rootsplit_count, parent_to_range);
    if (alpha > 0.) {
      // Last line of the section on EM in doc/tex.
      score_history[em_idx] += m_tilde_for_positive_alpha.dot(sbn_parameters);
    }
    // Return if we've converged according to score.
    if (em_idx > 0) {
      double scaled_score_improvement =
          (score_history[em_idx] - score_history[em_idx - 1]) /
          fabs(score_history[em_idx - 1]);
      // To monitor correctness of EM, we check to ensure that the score is
      // monotonically increasing (modulo numerical instability).
      // SHJ: -EPS is too small, I noticed the assertion failure for
      // scaled_score_improvement of -6e-16. Using ERR_TOLERANCE.
      Assert(scaled_score_improvement > -ERR_TOLERANCE, "Score function decreased.");
      if (fabs(scaled_score_improvement) < score_epsilon) {
        std::cout << "EM converged according to normalized score improvement < "
                  << score_epsilon << "." << std::endl;
        score_history.resize(em_idx + 1);
        break;
      }
    }
    ++progress_bar;
    progress_bar.display();
  }  // End of EM loop.
  progress_bar.done();
  NumericalUtils::ReportFloatingPointEnvironmentExceptions("|After EM|");
  return score_history;
}

bool SBNProbability::IsInSBNSupport(
    const RootedIndexerRepresentation& rooted_representation,
    size_t out_of_support_sentinel_value) {
  for (size_t idx : rooted_representation) {
    // Our convention is that out_of_support_sentinel_value is one more than the maximum
    // allowed PCSP index, so here we check the index is reasonable.
    Assert(idx <= out_of_support_sentinel_value,
           "Rooted tree index is greater than maximum permitted.");
    if (idx == out_of_support_sentinel_value) {
      return false;
    }
  }
  return true;
};

double SBNProbability::ProbabilityOfSingle(
    const EigenConstVectorXdRef sbn_parameters,
    const RootedIndexerRepresentation& rooted_representation) {
  size_t sbn_parameter_count = sbn_parameters.size();
  return IsInSBNSupport(rooted_representation, sbn_parameter_count)
             ? exp(SumOf(sbn_parameters, rooted_representation, 0.))
             : 0.;
}

double SBNProbability::ProbabilityOfSingle(
    const EigenConstVectorXdRef sbn_parameters,
    const UnrootedIndexerRepresentation& indexer_representation) {
  size_t sbn_parameter_count = sbn_parameters.size();
  double log_total_probability = DOUBLE_NEG_INF;
  for (const auto& rooted_representation : indexer_representation) {
    log_total_probability = NumericalUtils::LogAdd(
        log_total_probability,
        IsInSBNSupport(rooted_representation, sbn_parameter_count)
            ? SumOf(sbn_parameters, rooted_representation, 0.)
            : DOUBLE_NEG_INF);
  }
  return exp(log_total_probability);
}

// This code is duplicated for the unrooted case below.
EigenVectorXd SBNProbability::ProbabilityOfCollection(
    const EigenConstVectorXdRef sbn_parameters,
    const std::vector<RootedIndexerRepresentation>& indexer_representations) {
  std::function<double(const RootedIndexerRepresentation&)> f =
      [sbn_parameters](const RootedIndexerRepresentation& indexer_representation) {
        return ProbabilityOfSingle(sbn_parameters, indexer_representation);
      };
  return MapFromStdVectorToEigenVectorXd(indexer_representations, f);
}

// This code is duplicated for the rooted case above.
EigenVectorXd SBNProbability::ProbabilityOfCollection(
    const EigenConstVectorXdRef sbn_parameters,
    const std::vector<UnrootedIndexerRepresentation>& indexer_representations) {
  const size_t topology_count = indexer_representations.size();
  EigenVectorXd results(topology_count);
  for (size_t topology_idx = 0; topology_idx < topology_count; ++topology_idx) {
    results[topology_idx] =
        ProbabilityOfSingle(sbn_parameters, indexer_representations[topology_idx]);
  }
  return results;
}
