// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef GP_OPERATION_HPP_
#define GP_OPERATION_HPP_

#include <iostream>
#include <variant>

namespace GPOperations {

struct Zero {
  size_t idx;
};

struct WeightedSumAccumulate {
  size_t dest_idx;
  size_t weight_idx;
  size_t src_idx;
};

struct Multiply {
  size_t dest_idx;
  size_t src1_idx;
  size_t src2_idx;
};

struct Likelihood {
  size_t dest_idx;
  size_t src1_idx;
  size_t src2_idx;
};

struct EvolveRootward {
  size_t dest_idx;
  size_t src_idx;
  size_t branch_length_idx;
};

struct EvolveLeafward {
  size_t dest_idx;
  size_t src_idx;
  size_t branch_length_idx;
};

struct OptimizeRootward {
  size_t plv_dest_idx;
  size_t leafward_idx;
  size_t rootward_idx;
  size_t branch_length_idx;
};

struct OptimizeLeafward {
  size_t plv_dest_idx;
  size_t leafward_idx;
  size_t rootward_idx;
  size_t branch_length_idx;
};

struct UpdateSBNProbabilities {
  size_t start_idx;
  size_t stop_idx;
};
}  // namespace GPOperations

using GPOperation =
    std::variant<GPOperations::Zero, GPOperations::WeightedSumAccumulate>;
//                 GPOperations::Multiply, GPOperations::Likelihood,
//                 GPOperations::EvolveRootward, GPOperations::EvolveLeafward,
//                 GPOperations::OptimizeRootward, GPOperations::OptimizeLeafward,
//                 GPOperations::UpdateSBNProbabilities>;

struct GPOperationOstream {
  std::ostream& os_;

  GPOperationOstream(std::ostream& os) : os_{os} {}

  template <typename... Args>
  void FoldPrint(Args&&... args) {
    (os_ << ... << std::forward<Args>(args)) << '\n';
  }

  void operator()(const GPOperations::Zero& operation) {
    FoldPrint("Zero", operation.idx);
  }
  void operator()(const GPOperations::WeightedSumAccumulate& operation) {
    os_ << "WeightedSumAccumulate";
  }
};

std::ostream& operator<<(std::ostream& os, GPOperation const& v) {
  std::visit(GPOperationOstream{os}, v);
  return os;
}

/*
## Assumptions

* we have all PLVs in a data structure such that `plv[idx]` is the `idx`th PLV
* we have branch lengths in a vector `branch_lengths` indexed by PCSPs
* we have likelihoods in a vector `likelihoods` indexed by PCSPs
* we have a vector of SBN probabilities `q` indexed by PCSPs, such that the children of
a given parent are contiguous

### Desired operations

* `Zero(idx)` zeroes out the PLV at `idx`
* `WeightedSumAccumulate(dest_idx, weight_idx, src_idx)` does `plv[dest_idx] +=
q[weight_idx] * plv[src_idx]`
* `Multiply(dest_idx, src1_idx, src2_idx)` does `plv[dest_idx] = plv[src1_idx] circ
plv[src2_idx]`
* `Likelihood(dest_idx, src1_idx, src2_idx)` stores the "dot product" (flattened,
denoted `.`) of `plv[src1_idx]` and `plv[src2_idx]` in `likelihoods[dest_idx]`
* `EvolveRootward(dest_idx, src_idx, branch_length_idx)` does `plv[dest_idx] =
P(branch_lengths[branch_length_idx]) plv[src_idx]`
* `EvolveLeafward(dest_idx, src_idx, branch_length_idx)` does `plv[dest_idx] =
P'(branch_lengths[branch_length_idx]) plv[src_idx]`
* `OptimizeRootward(plv_dest_idx, leafward_idx, rootward_idx, branch_length_idx)` finds
the optimal `branch_length` for the likelihood `plv[rootward_idx] . P(branch_length)
plv[leafward_idx]`, starting optimization at `branch_lengths[branch_length_idx]`, and
stores the PLV for `P(branch_length) plv[leafward_idx]` in `plv[dest_idx]` for the
optimal branch length as well as the optimal branch length at
`branch_lengths[branch_length_idx]`
* `OptimizeLeafward(plv_dest_idx, leafward_idx, rootward_idx, branch_length_idx)` finds
the optimal `branch_length` for the likelihood `(P'(branch_length) plv[rootward_idx]) .
plv[leafward_idx]`, starting optimization at `branch_lengths[branch_length_idx]`, and
stores the PLV for `P'(branch_length) plv[rootward_idx]` in `plv[dest_idx]` for the
optimal branch length as well as the optimal branch length at
`branch_lengths[branch_length_idx]`
* `UpdateSBNProbabilities(start_idx, end_idx)` performs `eq:SBNUpdates`. That is, let
`total` be the sum of `likelihoods[idx]` for all `idx` in `start_idx <= idx < end_idx`.
Now let `q[idx] = likelihoods[idx]/total` for all `idx` in `start_idx <= idx < end_idx`.

Implement these as a
[std::variant](https://arne-mertz.de/2018/05/modern-c-features-stdvariant-and-stdvisit/)
, and do graph traversal to build a list of them.

*/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("TODO") {}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // GP_OPERATION_HPP_
