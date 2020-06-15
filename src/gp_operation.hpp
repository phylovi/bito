// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// These operations are just declarations. How we process them is a separate matter.

#ifndef GP_OPERATION_HPP_
#define GP_OPERATION_HPP_

#include <iostream>
#include <variant>
#include <vector>
#include "sugar.hpp"

using StringSizePairVector = std::vector<std::pair<std::string, size_t>>;

// Assume we have:
// * all PLVs in a data structure such that `plv[idx]` is the `idx`th PLV
// * branch lengths in a vector `branch_lengths` indexed by PCSPs
// * log likelihoods in a vector `log_likelihoods` indexed by PCSPs
// * a vector of SBN probabilities `q` indexed by PCSPs, such that the children of a
// given parent are contiguous
// * a collection of weights from site pattern compression

namespace GPOperations {

// We use the convention that `src_` and `dest_` indices always index PLVs.

// Zero out the PLV at `dest_idx`.
struct Zero {
  size_t dest_idx;
  StringSizePairVector guts() const { return {{"dest_idx", dest_idx}}; }
};

// Set the PLV at `dest_idx` to be the stationary distribution at every site.
struct SetToStationaryDistribution {
  size_t dest_idx;
  StringSizePairVector guts() const { return {{"dest_idx", dest_idx}}; }
};

// Perform `plv[dest_idx] += q[q_idx] * plv[src_idx]`
struct WeightedSumAccumulate {
  size_t dest_idx;
  size_t pcsp_idx;
  size_t src_idx;
  StringSizePairVector guts() const {
    return {{"dest_idx", dest_idx}, {"pcsp_idx", pcsp_idx}, {"src_idx", src_idx}};
  }
};

// Perform `plv[dest_idx] += q[q_idx] * plv[src_idx]`
struct MarginalLikelihood {
  size_t stationary_idx;
  size_t pcsp_idx;
  size_t p_idx;
  StringSizePairVector guts() const {
    return {{"r_idx", stationary_idx}, {"pcsp_idx", pcsp_idx}, {"p_idx", p_idx}};
  }
};

// Componentwise multiplication: `plv[dest_idx] = plv[src1_idx] o plv[src2_idx]`
struct Multiply {
  size_t dest_idx;
  size_t src1_idx;
  size_t src2_idx;
  StringSizePairVector guts() const {
    return {{"dest_idx", dest_idx}, {"src1_idx", src1_idx}, {"src2_idx", src2_idx}};
  }
};

// Stores the likelihood of `plv[child_idx]` and `plv[parent_idx]`, incorporating site
// pattern weights, in `log_likelihoods[dest_idx]`
struct Likelihood {
  size_t dest_idx;
  size_t child_idx;
  size_t parent_idx;
  StringSizePairVector guts() const {
    return {{"dest_idx", dest_idx}, {"child_idx", child_idx}, {"parent_idx", parent_idx}};
  }
};

// Perform `plv[dest_idx] = P(branch_lengths[branch_length_idx]) plv[src_idx]`
struct EvolveRootward {
  size_t dest_idx;
  size_t src_idx;
  size_t branch_length_idx;
  StringSizePairVector guts() const {
    return {{"dest_idx", dest_idx},
            {"src_idx", src_idx},
            {"branch_length_idx", branch_length_idx}};
  }
};

// Perform `plv[dest_idx] = P'(branch_lengths[branch_length_idx]) plv[src_idx]`
struct EvolveLeafward {
  size_t dest_idx;
  size_t src_idx;
  size_t branch_length_idx;
  StringSizePairVector guts() const {
    return {{"dest_idx", dest_idx},
            {"src_idx", src_idx},
            {"branch_length_idx", branch_length_idx}};
  }
};

// Finds the optimal `branch_length` for the likelihood of
// `plv[rootward_idx]` and `P(branch_length) plv[leafward_idx]`,
// * starting optimization at `branch_lengths[branch_length_idx]`,
// * storing the PLV for `P(branch_length) plv[leafward_idx]` in `plv[dest_idx]` for the
// optimal branch length
// * storing log likelihood at `log_likelihoods[branch_length_idx]`
// * storing optimal branch length at `branch_lengths[branch_length_idx]`
struct OptimizeBranchLength {
  size_t dest_idx;
  size_t leafward_idx;
  size_t rootward_idx;
  size_t branch_length_idx;
  StringSizePairVector guts() const {
    return {{"dest_idx", dest_idx},
            {"leafward_idx", leafward_idx},
            {"rootward_idx", rootward_idx},
            {"branch_length_idx", branch_length_idx}};
  }
};

// Performs `eq:SBNUpdates`. That is, let `total` be the sum of `log_likelihoods[idx]`
// for all `idx` in `start_idx <= idx < stop_idx`. Now let `q[idx] =
// log_likelihoods[idx]/total` for all `idx` in `start_idx <= idx < stop_idx`.
struct UpdateSBNProbabilities {
  size_t start_idx;
  size_t stop_idx;
  StringSizePairVector guts() const {
    return {{"start_idx", start_idx}, {"stop_idx", stop_idx}};
  }
};
}  // namespace GPOperations

using GPOperation =
    std::variant<GPOperations::Zero, GPOperations::SetToStationaryDistribution,
                 GPOperations::WeightedSumAccumulate, GPOperations::Multiply,
                 GPOperations::Likelihood, GPOperations::EvolveRootward,
                 GPOperations::EvolveLeafward, GPOperations::OptimizeBranchLength,
                 GPOperations::UpdateSBNProbabilities,
                 GPOperations::MarginalLikelihood>;

using GPOperationVector = std::vector<GPOperation>;

struct GPOperationOstream {
  std::ostream& os_;

  GPOperationOstream(std::ostream& os) : os_{os} {}

  void operator()(const GPOperations::Zero& operation) {
    os_ << "Zero" << operation.guts();
  }
  void operator()(const GPOperations::SetToStationaryDistribution& operation) {
    os_ << "SetToStationaryDistribution" << operation.guts();
  }
  void operator()(const GPOperations::WeightedSumAccumulate& operation) {
    os_ << "WeightedSumAccumulate" << operation.guts();
  }
  void operator()(const GPOperations::MarginalLikelihood& operation) {
    os_ << "MarginalLikelihood" << operation.guts();
  }
  void operator()(const GPOperations::Multiply& operation) {
    os_ << "Multiply" << operation.guts();
  }
  void operator()(const GPOperations::Likelihood& operation) {
    os_ << "Likelihood" << operation.guts();
  }
  void operator()(const GPOperations::EvolveRootward& operation) {
    os_ << "EvolveRootward" << operation.guts();
  }
  void operator()(const GPOperations::EvolveLeafward& operation) {
    os_ << "EvolveLeafward" << operation.guts();
  }
  void operator()(const GPOperations::OptimizeBranchLength& operation) {
    os_ << "OptimizeBranchLength" << operation.guts();
  }
  void operator()(const GPOperations::UpdateSBNProbabilities& operation) {
    os_ << "UpdateSBNProbabilities" << operation.guts();
  }
};

std::ostream& operator<<(std::ostream& os, GPOperation const& operation);
std::ostream& operator<<(std::ostream& os, GPOperationVector const& operation_vector);

#endif  // GP_OPERATION_HPP_
