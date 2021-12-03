// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "topology_sampler.hpp"
#include "numerical_utils.hpp"

Node::NodePtr TopologySampler::SampleTopology(const Input &input) const {
    // Start by sampling a rootsplit.
    size_t rootsplit_index =
        SampleIndex(input, std::pair<size_t, size_t>(0, input.RootsplitCount()));
    const Bitset &rootsplit = input.RootsplitsAt(rootsplit_index);
    auto topology =
        input.IsRooted() ? SampleTopology(input, rootsplit) : SampleTopology(input, rootsplit)->Deroot();
    topology->Polish();
    return topology;
}

// Sample an integer index in [range.first, range.second) according to
// SBNParameters.
size_t TopologySampler::SampleIndex(const Input &input, Range range) const {
  const auto &[start, end] = range;
  Assert(start < end && static_cast<Eigen::Index>(end) <= input.SBNParameters().size(),
          "SampleIndex given an invalid range.");
  // We do not want to overwrite sbn_parameters so we make a copy.
  EigenVectorXd sbn_parameters_subrange = input.SBNParameters().segment(start, end - start);
  NumericalUtils::ProbabilityNormalizeInLog(sbn_parameters_subrange);
  NumericalUtils::Exponentiate(sbn_parameters_subrange);
  std::discrete_distribution<> distribution(sbn_parameters_subrange.begin(),
                                            sbn_parameters_subrange.end());
  // We have to add on range.first because we have taken a slice of the full
  // array, and the sampler treats the beginning of this slice as zero.
  auto result =
      start + static_cast<size_t>(distribution(mersenne_twister_.GetGenerator()));
  Assert(result < end, "SampleIndex sampled a value out of range.");
  return result;
}

// The input to this function is a parent subsplit (of length 2n).
Node::NodePtr TopologySampler::SampleTopology(const Input &input, const Bitset &parent_subsplit) const {
  auto process_subsplit = [this, &input](const Bitset &parent) {
    auto singleton_option = parent.SubsplitGetClade(1).SingletonOption();
    if (singleton_option) {
      return Node::Leaf(*singleton_option);
    }  // else
    auto child_index = SampleIndex(input, input.ParentToRangeAt(parent));
    return SampleTopology(input, input.IndexToChildAt(child_index));
  };
  return Node::Join(process_subsplit(parent_subsplit),
                    process_subsplit(parent_subsplit.SubsplitRotate()));
}