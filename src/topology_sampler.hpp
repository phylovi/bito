// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "node.hpp"
#include "mersenne_twister.hpp"
#include "eigen_sugar.hpp"

class TopologySampler {
public:

  class Input {
  public:
    virtual bool IsRooted() const = 0;
    virtual EigenConstVectorXdRef SBNParameters() const = 0;
    virtual size_t RootsplitCount() const = 0;
    virtual const Bitset &RootsplitsAt(size_t rootsplit_idx) const = 0;
    virtual const SizePair &ParentToRangeAt(const Bitset &parent) const = 0;
    virtual const Bitset &IndexToChildAt(size_t child_idx) const = 0;
  };

  Node::NodePtr SampleTopology(const Input& input) const;
  inline void SetSeed(uint64_t seed) { mersenne_twister_.SetSeed(seed); }

private:
  using Range = std::pair<size_t, size_t>;

  size_t SampleIndex(const Input &input, Range range) const;
  Node::NodePtr SampleTopology(const Input &input, const Bitset &parent_subsplit) const;
  
  MersenneTwister mersenne_twister_;
};