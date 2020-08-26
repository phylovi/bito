// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// These operations are just declarations. We process them with the GPEngine.

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

// Zero out the PLV at `dest_`.
struct Zero {
  constexpr explicit Zero(size_t dest) : dest_{dest} {}
  size_t dest_;
  StringSizePairVector guts() const { return {{"dest_", dest_}}; }
};

// Set the PLV at `dest_` to be the stationary distribution at every site.
struct SetToStationaryDistribution {
  constexpr explicit SetToStationaryDistribution(size_t dest) : dest_{dest} {}
  size_t dest_;
  StringSizePairVector guts() const { return {{"dest_", dest_}}; }
};

// Set transition_matrix_ using branch_length(gpcsp_) then,
// perform `plv[dest_] += q[gpcsp_] * transition_matrix_ * plv[src_]`
struct IncrementWithWeightedEvolvedPLV {
  constexpr IncrementWithWeightedEvolvedPLV(size_t dest, size_t gpcsp, size_t src)
      : dest_{dest}, gpcsp_{gpcsp}, src_{src} {}
  size_t dest_;
  size_t gpcsp_;
  size_t src_;
  StringSizePairVector guts() const {
    return {{"dest_", dest_}, {"gpcsp_", gpcsp_}, {"src_", src_}};
  }
};

// Increment log marginal likelihood with the log likelihood at rootsplit rootsplit_
// with a leafward PLV p_ and a stationary distribution at stationary_.
struct IncrementMarginalLikelihood {
  constexpr IncrementMarginalLikelihood(size_t stationary, size_t rootsplit, size_t p)
      : stationary_{stationary}, rootsplit_{rootsplit}, p_{p} {}
  size_t stationary_;
  size_t rootsplit_;
  size_t p_;
  StringSizePairVector guts() const {
    return {{"stationary_", stationary_}, {"rootsplit_", rootsplit_}, {"p_", p_}};
  }
};

// Componentwise multiplication: `plv[dest_] = plv[src1_] o plv[src2_]`
struct Multiply {
  constexpr Multiply(size_t dest, size_t src1, size_t src2)
      : dest_{dest}, src1_{src1}, src2_{src2} {}
  size_t dest_;
  size_t src1_;
  size_t src2_;
  StringSizePairVector guts() const {
    return {{"dest_", dest_}, {"src1_", src1_}, {"src2_", src2_}};
  }
};

// Stores the likelihood of `plv[child_]` and `plv[parent_]` with branch length
// branch_lengths[dest_], incorporating site pattern weights, in
// `log_likelihoods[dest_]`
struct Likelihood {
  constexpr Likelihood(size_t dest, size_t child, size_t parent)
      : dest_{dest}, child_{child}, parent_{parent} {}
  size_t dest_;
  size_t child_;
  size_t parent_;
  StringSizePairVector guts() const {
    return {{"dest_", dest_}, {"child_", child_}, {"parent_", parent_}};
  }
};

// Finds the optimal `branch_length` for the likelihood of
// `plv[rootward_]` and `P(branch_length) plv[leafward_]`,
// * starting optimization at `branch_lengths[branch_length_]`,
// * storing the PLV for `P(branch_length) plv[leafward_]` in `plv[dest_]` for the
// optimal branch length
// * storing log likelihood at `log_likelihoods[branch_length_]`
// * storing optimal branch length at `branch_lengths[branch_length_]`
struct OptimizeBranchLength {
  constexpr OptimizeBranchLength(size_t leafward, size_t rootward, size_t gpcsp)
      : leafward_{leafward}, rootward_{rootward}, gpcsp_{gpcsp} {}
  size_t leafward_;
  size_t rootward_;
  size_t gpcsp_;
  StringSizePairVector guts() const {
    return {{"leafward_", leafward_}, {"rootward_", rootward_}, {"gpcsp_", gpcsp_}};
  }
};

// Assumption: log_likelihoods_ have been updated on [op.start_, op.stop_).
// Performs `eq:SBNUpdates`. That is, let `total` be the log sum of
// `log_likelihoods[idx]` for all `idx` in `start_ <= idx < stop_`. Now let
// `q[idx] = exp(log_likelihoods[idx] - total)` for all `idx` in `start_ <= idx <
// stop_`.
// Note that this operation modifies our log_likelihoods in place by normalizing them
// across children of a parent. Thus they are no longer valid.
struct UpdateSBNProbabilities {
  constexpr UpdateSBNProbabilities(size_t start, size_t stop)
      : start_{start}, stop_{stop} {}
  size_t start_;
  size_t stop_;
  StringSizePairVector guts() const { return {{"start_", start_}, {"stop_", stop_}}; }
};

// This operation sets the rescaling amount for the PLV in op.dest_ to be the minimum of
// that for all of the PLVs in op.src_vector_. We do this so that we can sum over
// partial likelihood vectors after rescaling each one so that it is on the same scale
// as dest. Note that we want the minimum here because we want to preserve accuracy for
// the PLV with the highest likelihood (which corresponds to the least amount of
// rescaling).
struct PrepForMarginalization {
  PrepForMarginalization(size_t dest, SizeVector src_vector)
      : dest_{dest}, src_vector_{std::move(src_vector)} {}
  size_t dest_;
  SizeVector src_vector_;
  std::pair<std::pair<std::string, size_t>, std::pair<std::string, SizeVector>> guts()
      const {
    return {{"dest_", dest_}, {"src_vector_", src_vector_}};
  }
};

}  // namespace GPOperations

using GPOperation =
    std::variant<GPOperations::Zero, GPOperations::SetToStationaryDistribution,
                 GPOperations::IncrementWithWeightedEvolvedPLV, GPOperations::Multiply,
                 GPOperations::Likelihood, GPOperations::OptimizeBranchLength,
                 GPOperations::UpdateSBNProbabilities,
                 GPOperations::IncrementMarginalLikelihood,
                 GPOperations::PrepForMarginalization>;

using GPOperationVector = std::vector<GPOperation>;

// The purpose of this visitor class is to accumulate the
// things-that-need-preparation-for-marginalization and build them into a
// PrepForMarginalization (see PrepForMarginalizationOfOperations implementation).
struct PrepForMarginalizationVisitor {
  std::optional<size_t> dest_ = std::nullopt;
  SizeVector src_vector;

  explicit PrepForMarginalizationVisitor(const GPOperationVector& operations) {
    for (const auto& operation : operations) {
      std::visit(*this, operation);
    }
  }

  void operator()(const GPOperations::IncrementWithWeightedEvolvedPLV& op) {
    if (dest_) {
      Assert(*dest_ == op.dest_, "dest_ mismatch in PrepForMarginalizationVisitor");
    } else {
      dest_ = op.dest_;
    }
    src_vector.push_back(op.src_);
  }
  // Do nothing for the rest of the operations.
  void operator()(const GPOperations::Zero&) {}                         // NOLINT
  void operator()(const GPOperations::SetToStationaryDistribution&) {}  // NOLINT
  void operator()(const GPOperations::IncrementMarginalLikelihood&) {}  // NOLINT
  void operator()(const GPOperations::Multiply&) {}                     // NOLINT
  void operator()(const GPOperations::Likelihood&) {}                   // NOLINT
  void operator()(const GPOperations::OptimizeBranchLength&) {}         // NOLINT
  void operator()(const GPOperations::UpdateSBNProbabilities&) {}       // NOLINT
  void operator()(const GPOperations::PrepForMarginalization&) {}       // NOLINT

  GPOperations::PrepForMarginalization ToPrepForMarginalization() {
    Assert(dest_, "Nothing to prep in ToPrepForMarginalization");
    return GPOperations::PrepForMarginalization{*dest_, src_vector};
  }
};

namespace GPOperations {
PrepForMarginalization PrepForMarginalizationOfOperations(
    const GPOperationVector& operations);
};  // namespace GPOperations

struct GPOperationOstream {
  std::ostream& os_;

  explicit GPOperationOstream(std::ostream& os) : os_{os} {}

  void operator()(const GPOperations::Zero& operation) {
    os_ << "Zero" << operation.guts();
  }
  void operator()(const GPOperations::SetToStationaryDistribution& operation) {
    os_ << "SetToStationaryDistribution" << operation.guts();
  }
  void operator()(const GPOperations::IncrementWithWeightedEvolvedPLV& operation) {
    os_ << "IncrementWithWeightedEvolvedPLV" << operation.guts();
  }
  void operator()(const GPOperations::IncrementMarginalLikelihood& operation) {
    os_ << "IncrementMarginalLikelihood" << operation.guts();
  }
  void operator()(const GPOperations::Multiply& operation) {
    os_ << "Multiply" << operation.guts();
  }
  void operator()(const GPOperations::Likelihood& operation) {
    os_ << "Likelihood" << operation.guts();
  }
  void operator()(const GPOperations::OptimizeBranchLength& operation) {
    os_ << "OptimizeBranchLength" << operation.guts();
  }
  void operator()(const GPOperations::UpdateSBNProbabilities& operation) {
    os_ << "UpdateSBNProbabilities" << operation.guts();
  }
  void operator()(const GPOperations::PrepForMarginalization& operation) {
    os_ << "PrepForMarginalization" << operation.guts();
  }
};

std::ostream& operator<<(std::ostream& os, GPOperation const& operation);
std::ostream& operator<<(std::ostream& os, GPOperationVector const& operation_vector);

#endif  // GP_OPERATION_HPP_
