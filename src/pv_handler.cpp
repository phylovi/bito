// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "pv_handler.hpp"

// ** Resize

template <class PVType, class PVTypeEnum>
void PVNodeHandler<PVType, PVTypeEnum>::Resize(const size_t new_node_count,
                                               const size_t new_node_alloc) {
  const size_t old_pv_count = GetPVCount();
  SetNodeCount(new_node_count);
  SetAllocatedNodeCount(new_node_alloc);
  // Allocate mmapped data block.
  mmapped_master_pvs_.Resize(GetAllocatedPVCount() * pattern_count_);
  // Subdivide mmapped data in individual PVs.
  pvs_ = mmapped_master_pvs_.Subdivide(GetAllocatedPVCount());
  // Initialize new work space.
  Assert((pvs_.back().rows() == MmappedNucleotidePLV::base_count_) &&
             (pvs_.back().cols() == static_cast<Eigen::Index>(pattern_count_)) &&
             (size_t(pvs_.size()) == GetAllocatedPVCount()),
         "Didn't get the right shape of PVs out of Subdivide.");
  for (size_t i = old_pv_count; i < GetPaddedPVCount(); i++) {
    pvs_.at(i).setZero();
  }
}

template <class PVType, class PVTypeEnum>
void PVEdgeHandler<PVType, PVTypeEnum>::Resize(const size_t new_edge_count,
                                               const size_t new_edge_alloc) {
  const size_t old_pv_count = GetPVCount();
  SetEdgeCount(new_edge_count);
  SetAllocatedEdgeCount(new_edge_alloc);
  // Allocate mmapped data block.
  mmapped_master_pvs_.Resize(GetAllocatedPVCount() * pattern_count_);
  // Subdivide mmapped data in individual PVs.
  pvs_ = mmapped_master_pvs_.Subdivide(GetAllocatedPVCount());
  // Initialize new work space.
  Assert((pvs_.back().rows() == MmappedNucleotidePLV::base_count_) &&
             (pvs_.back().cols() == static_cast<Eigen::Index>(pattern_count_)) &&
             (size_t(pvs_.size()) == GetAllocatedPVCount()),
         "Didn't get the right shape of PVs out of Subdivide.");
  for (size_t i = old_pv_count; i < GetPaddedPVCount(); i++) {
    pvs_.at(i).setZero();
  }
}

template <class PVType, class PVTypeEnum>
void PartialVectorHandler<PVType, PVTypeEnum>::Reindex(const Reindexer pv_reindexer) {
  Reindexer::ReindexInPlace(pvs_, pv_reindexer, GetPVCount(), GetPV(GetPVCount()),
                            GetPV(GetPVCount() + 1));
}

template <class PVType, class PVTypeEnum>
Reindexer PartialVectorHandler<PVType, PVTypeEnum>::BuildPVReindexer(
    const Reindexer& node_reindexer, const size_t old_node_count,
    const size_t new_node_count) {
  node_count_ = new_node_count;
  Reindexer pv_reindexer(new_node_count * pv_count_per_node_);
  size_t new_pvs_idx = old_node_count * pv_count_per_node_;
  for (size_t i = 0; i < new_node_count; i++) {
    const NodeId old_node_idx = NodeId(i);
    const NodeId new_node_idx =
        NodeId(node_reindexer.GetNewIndexByOldIndex(old_node_idx.value_));
    for (const auto pv_type : typename PVTypeEnum::Iterator()) {
      // Either get input pv_index from old pvs, or get new pv_index (new data is
      // irrelevant, so just get next available index).
      size_t old_pv_idx;
      if (old_node_idx < old_node_count) {
        old_pv_idx = GetPVIndex(pv_type, old_node_idx, old_node_count);
      } else {
        old_pv_idx = new_pvs_idx;
        new_pvs_idx++;
      }
      const size_t new_pv_idx = GetPVIndex(pv_type, new_node_idx, new_node_count);
      pv_reindexer.SetReindex(old_pv_idx, new_pv_idx);
    }
  }
  Assert(pv_reindexer.IsValid(GetPVCount()), "PV Reindexer is not valid.");
  return pv_reindexer;
}

// ** Explicit Instantiation
template class PartialVectorHandler<PartialVectorType::PLVType,
                                    PartialVectorType::PLVTypeEnum>;
template class PartialVectorHandler<PartialVectorType::PSVType,
                                    PartialVectorType::PSVTypeEnum>;
