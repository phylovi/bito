// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// NNI Engine Key Index
// Map for storing important key indices for NNI Likelihood computation of proposed
// NNIs. Creates a mapping of specified NodeId, EdgeId, and PVIds from pre-NNI to
// post-NNI.

#pragma once

#include "wrapper_sugar.hpp"
#include "subsplit_dag_storage.hpp"

// ** Key Indexing
// These function finds mapping from current NNIs contained in DAG to proposed NNI.
enum class NNIEngineKeyIndex : size_t {
  Parent_Id,
  Child_Id,
  Edge,
  Parent_RHat,
  Parent_RFocal,
  Parent_PHatSister,
  Child_P,
  Child_PHatLeft,
  Child_PHatRight,
};
static const size_t NNIEngineKeyIndexCount = 9;
class NNIEngineKeyIndexEnum
    : public EnumWrapper<NNIEngineKeyIndex, size_t, NNIEngineKeyIndexCount,
                         NNIEngineKeyIndex::Parent_Id,
                         NNIEngineKeyIndex::Child_PHatRight> {};

using NNIEngineKeyIndexMap = NNIEngineKeyIndexEnum::Array<size_t>;
using NNIEngineKeyIndexMapPair = std::pair<NNIEngineKeyIndexMap, NNIEngineKeyIndexMap>;
using NNIEngineKeyIndexPairArray =
    std::array<std::pair<NNIEngineKeyIndex, NNIEngineKeyIndex>, 3>;

enum class NNIEngineNodeIndex : size_t {
  Grandparent,
  Parent,
  Child,
  Sister,
  LeftChild,
  RightChild,
};
static const size_t NNIEngineNodeIndexCount = 6;
class NNIEngineNodeIndexEnum
    : public EnumWrapper<NNIEngineNodeIndex, size_t, NNIEngineNodeIndexCount,
                         NNIEngineNodeIndex::Grandparent,
                         NNIEngineNodeIndex::RightChild> {};
using NNIEngineNodeIndexMap = NNIEngineNodeIndexEnum::Array<NodeId>;

enum class NNIEngineEdgeIndex : size_t {
  Parent,
  Central,
  Sister,
  LeftChild,
  RightChild
};
static const size_t NNIEngineEdgeIndexCount = 5;
class NNIEngineEdgeIndexEnum
    : public EnumWrapper<NNIEngineEdgeIndex, size_t, NNIEngineEdgeIndexCount,
                         NNIEngineEdgeIndex::Parent, NNIEngineEdgeIndex::RightChild> {};
using NNIEngineEdgeIndexMap = NNIEngineEdgeIndexEnum::Array<EdgeId>;
