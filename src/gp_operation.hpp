// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// These operations are just declarations. We process them with the GPEngine.

#pragma once

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

// ZeroPLV out the PLV at `dest_`.
struct ZeroPLV {
  constexpr explicit ZeroPLV(size_t dest) : dest_{dest} {}
  size_t dest_;
  StringSizePairVector guts() const { return {{"dest_", dest_}}; }
};

// Set the PLV at `dest_` to be the stationary distribution multiplied by
// rootsplit probability indexed by `root_gpcsp_idx` at every site.
struct SetToStationaryDistribution {
  constexpr explicit SetToStationaryDistribution(size_t dest, size_t root_gpcsp_idx)
      : dest_{dest}, root_gpcsp_idx_(root_gpcsp_idx) {}
  size_t dest_;
  size_t root_gpcsp_idx_;
  StringSizePairVector guts() const {
    return {{"dest_", dest_}, {"root_gpcsp_idx_", root_gpcsp_idx_}};
  }
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

// Reset the marginal likelihood before incrementing.
struct ResetMarginalLikelihood {
  constexpr ResetMarginalLikelihood() {}
  StringSizePairVector guts() const { return {}; }
};

// Increment log marginal likelihood with the log likelihood at rootsplit rootsplit_
// with a leafward PLV p_ and a stationary distribution multiplied by the rootsplit
// prior at stationary_times_prior_.
// #288 note that this also sets the per-rootsplit conditional likelihood. It deserves a
// better name that has something to do with rootsplits. We should also consider just
// rejiggering things to use the stationary distribution directly.
struct IncrementMarginalLikelihood {
  constexpr IncrementMarginalLikelihood(size_t stationary_times_prior, size_t rootsplit,
                                        size_t p)
      : stationary_times_prior_{stationary_times_prior}, rootsplit_{rootsplit}, p_{p} {}
  size_t stationary_times_prior_;
  size_t rootsplit_;
  size_t p_;
  StringSizePairVector guts() const {
    return {{"stationary_times_prior_", stationary_times_prior_},
            {"rootsplit_", rootsplit_},
            {"p_", p_}};
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

// #288 this deserves a better description, and perhaps a better name.
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
// starting optimization at `branch_lengths[branch_length_]`, and
// storing optimal branch length at `branch_lengths[branch_length_]`.
// #288 are we happy with definition of rootward and leafward?
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

using GPOperation = std::variant<
    GPOperations::ZeroPLV, GPOperations::SetToStationaryDistribution,
    GPOperations::IncrementWithWeightedEvolvedPLV, GPOperations::Multiply,
    GPOperations::Likelihood, GPOperations::OptimizeBranchLength,
    GPOperations::UpdateSBNProbabilities, GPOperations::ResetMarginalLikelihood,
    GPOperations::IncrementMarginalLikelihood, GPOperations::PrepForMarginalization>;

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
  void operator()(const GPOperations::ZeroPLV&) {}                      // NOLINT
  void operator()(const GPOperations::SetToStationaryDistribution&) {}  // NOLINT
  void operator()(const GPOperations::ResetMarginalLikelihood&) {}      // NOLINT
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
void AppendGPOperations(GPOperationVector& operations,
                        GPOperationVector&& new_operations);

PrepForMarginalization PrepForMarginalizationOfOperations(
    const GPOperationVector& operations);
};  // namespace GPOperations

struct GPOperationOstream {
  std::ostream& os_;

  explicit GPOperationOstream(std::ostream& os) : os_{os} {}

  void operator()(const GPOperations::ZeroPLV& operation) {
    os_ << "ZeroPLV" << operation.guts();
  }
  void operator()(const GPOperations::SetToStationaryDistribution& operation) {
    os_ << "SetToStationaryDistribution" << operation.guts();
  }
  void operator()(const GPOperations::IncrementWithWeightedEvolvedPLV& operation) {
    os_ << "IncrementWithWeightedEvolvedPLV" << operation.guts();
  }
  // #288 should this be TotalMarginalLikelihood or something?
  void operator()(const GPOperations::ResetMarginalLikelihood& operation) {
    os_ << "ResetMarginalLikelihood" << operation.guts();
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

