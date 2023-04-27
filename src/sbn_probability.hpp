// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Perform training of an SBN based on a sample of trees.
//
// We assume that readers are familiar with how the sbn_parameters_ vector is laid out:
// first probabilities of rootsplits, then conditional probabilities of PCSPs.

#pragma once

#include "eigen_sugar.hpp"
#include "sbn_maps.hpp"

namespace SBNProbability {

// The "SBN-SA" estimator described in the "Maximum Lower Bound Estimates" section of
// the 2018 NeurIPS paper.
void SimpleAverage(
    EigenVectorXdRef sbn_parameters,
    const UnrootedIndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range);

void SimpleAverage(
    EigenVectorXdRef sbn_parameters,
    const RootedIndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range);

// The "SBN-EM" estimator described in the "Expectation Maximization" section of
// the 2018 NeurIPS paper. Returns the sequence of scores (defined in the paper)
// obtained by the EM iterations.
EigenVectorXd ExpectationMaximization(
    EigenVectorXdRef sbn_parameters,
    const UnrootedIndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range, double alpha,
    size_t max_iter, double score_epsilon);

// Calculate the probability of an indexer_representation of a rooted topology.
double ProbabilityOfSingle(EigenConstVectorXdRef sbn_parameters,
                           const RootedIndexerRepresentation& indexer_representation);

// Calculate the probability of an indexer_representation of an unrooted topology.
double ProbabilityOfSingle(EigenConstVectorXdRef sbn_parameters,
                           const UnrootedIndexerRepresentation& indexer_representation);

// Calculate the probabilities of a collection of rooted indexer_representations.
EigenVectorXd ProbabilityOfCollection(
    EigenConstVectorXdRef sbn_parameters,
    const std::vector<RootedIndexerRepresentation>& indexer_representations);

// Calculate the probabilities of a collection of unrooted indexer_representations.
EigenVectorXd ProbabilityOfCollection(
    EigenConstVectorXdRef sbn_parameters,
    const std::vector<UnrootedIndexerRepresentation>& indexer_representations);

// This function performs in-place normalization of vec given by range when its values
// are in log space.
void ProbabilityNormalizeRangeInLog(EigenVectorXdRef vec,
                                    std::pair<size_t, size_t> range);
// Perform in-place normalization of vec when its values are in log space.
// We assume that vec is laid out like sbn_parameters (see top).
void ProbabilityNormalizeParamsInLog(EigenVectorXdRef vec, size_t rootsplit_count,
                                     const BitsetSizePairMap& parent_to_range);
bool IsInSBNSupport(const SizeVector& rooted_representation,
                    size_t out_of_support_sentinel_value);

// Take the sum of the entries of vec in indices plus starting_value.
double SumOf(EigenConstVectorXdRef vec, const SizeVector& indices,
             double starting_value);

}  // namespace SBNProbability
