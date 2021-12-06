// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "topology_sampler.hpp"
#include "numerical_utils.hpp"
#include "unrooted_sbn_instance.hpp"
#include "rooted_sbn_instance.hpp"
#include "subsplit_dag.hpp"

Node::NodePtr TopologySampler::SampleTopology(const Input &input, bool is_rooted) const {
    // Start by sampling a rootsplit.
    size_t rootsplit_index =
        SampleIndex(input, std::pair<size_t, size_t>(0, input.RootsplitCount()));
    const Bitset &rootsplit = input.RootsplitsAt(rootsplit_index);
    auto topology =
        is_rooted ? SampleTopology(input, rootsplit) : SampleTopology(input, rootsplit)->Deroot();
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
  auto singleton_option = parent_subsplit.SubsplitGetClade(0).SingletonOption();
  if (singleton_option && parent_subsplit.SubsplitGetClade(1).None()) {
    return Node::Leaf(*singleton_option);
  }
  singleton_option = parent_subsplit.SubsplitGetClade(1).SingletonOption();
  if (singleton_option && parent_subsplit.SubsplitGetClade(0).None()) {
    return Node::Leaf(*singleton_option);
  }
  auto process_subsplit = [this, &input](const Bitset &parent) {
    auto child_index = SampleIndex(input, input.ParentToRangeAt(parent));
    return SampleTopology(input, input.IndexToChildAt(child_index));
  };
  return Node::Join(process_subsplit(parent_subsplit),
                    process_subsplit(parent_subsplit.SubsplitRotate()));
}

size_t UnrootedSBNInstanceSamplerInput::RootsplitCount() const {
  return inst_.SBNSupport().RootsplitCount();
}

EigenConstVectorXdRef UnrootedSBNInstanceSamplerInput::SBNParameters() const {
  return inst_.SBNParameters();
}

const Bitset& UnrootedSBNInstanceSamplerInput::RootsplitsAt(size_t rootsplit_idx) const {
  return inst_.SBNSupport().RootsplitsAt(rootsplit_idx);
}

const SizePair& UnrootedSBNInstanceSamplerInput::ParentToRangeAt(const Bitset &parent) const {
  return inst_.SBNSupport().ParentToRangeAt(parent);
}

const Bitset& UnrootedSBNInstanceSamplerInput::IndexToChildAt(size_t child_idx) const {
  return inst_.SBNSupport().IndexToChildAt(child_idx);
}

size_t RootedSBNInstanceSamplerInput::RootsplitCount() const {
  return inst_.SBNSupport().RootsplitCount();
}

EigenConstVectorXdRef RootedSBNInstanceSamplerInput::SBNParameters() const {
  return inst_.SBNParameters();
}

const Bitset& RootedSBNInstanceSamplerInput::RootsplitsAt(size_t rootsplit_idx) const {
  return inst_.SBNSupport().RootsplitsAt(rootsplit_idx);
}

const SizePair& RootedSBNInstanceSamplerInput::ParentToRangeAt(const Bitset &parent) const {
  return inst_.SBNSupport().ParentToRangeAt(parent);
}

const Bitset& RootedSBNInstanceSamplerInput::IndexToChildAt(size_t child_idx) const {
  return inst_.SBNSupport().IndexToChildAt(child_idx);
}

size_t SubsplitDAGSamplerInput::RootsplitCount() const {
  return dag_.RootsplitCount();
}

EigenConstVectorXdRef SubsplitDAGSamplerInput::SBNParameters() const {
  return sbn_parameters_;
}

const Bitset& SubsplitDAGSamplerInput::RootsplitsAt(size_t rootsplit_idx) const {
  return dag_.GetDAGNode(dag_.RootsplitIds().at(rootsplit_idx))->GetBitset();;
}

const SizePair& SubsplitDAGSamplerInput::ParentToRangeAt(const Bitset &parent) const {
  return dag_.ParentToRange().at(parent);
}

const Bitset& SubsplitDAGSamplerInput::IndexToChildAt(size_t child_idx) const {
  for (auto&& i : dag_.DAGEdges()) {
    if (i.second == child_idx) return dag_.GetDAGNode(i.first.second)->GetBitset();
  }
  Failwith("Edge not found");
}
