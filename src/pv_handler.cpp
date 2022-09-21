// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "pv_handler.hpp"

// ** Resize

template <class PVType, class PVTypeEnum, class DAGElementId>
void PartialVectorHandler<PVType, PVTypeEnum, DAGElementId>::Resize(
    const size_t new_element_count, const size_t new_element_alloc) {
  const size_t old_pv_count = GetPVCount();
  element_count_ = new_element_count;
  element_alloc_ = new_element_alloc;
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

template <class PVType, class PVTypeEnum, class DAGElementId>
void PartialVectorHandler<PVType, PVTypeEnum, DAGElementId>::Reindex(
    const Reindexer pv_reindexer) {
  Reindexer::ReindexInPlace(pvs_, pv_reindexer, GetPVCount(), GetPV(GetPVCount()),
                            GetPV(GetPVCount() + 1));
}

template <class PVType, class PVTypeEnum, class DAGElementId>
Reindexer PartialVectorHandler<PVType, PVTypeEnum, DAGElementId>::BuildPVReindexer(
    const Reindexer& element_reindexer, const size_t old_element_count,
    const size_t new_element_count) {
  element_count_ = new_element_count;
  Reindexer pv_reindexer(new_element_count * pv_count_per_element_);
  PVId new_pvs_idx = PVId(old_element_count * pv_count_per_element_);
  for (size_t i = 0; i < new_element_count; i++) {
    const DAGElementId old_element_idx = DAGElementId(i);
    const DAGElementId new_element_idx =
        DAGElementId(element_reindexer.GetNewIndexByOldIndex(old_element_idx.value_));
    for (const auto pv_type : typename PVTypeEnum::Iterator()) {
      // Either get input pv_index from old pvs, or get new pv_index (new data is
      // irrelevant, so just get next available index).
      PVId old_pv_idx;
      if (old_element_idx < old_element_count) {
        old_pv_idx = GetPVIndex(pv_type, old_element_idx, old_element_count);
      } else {
        old_pv_idx = new_pvs_idx;
        new_pvs_idx++;
      }
      const PVId new_pv_idx = GetPVIndex(pv_type, new_element_idx, new_element_count);
      pv_reindexer.SetReindex(old_pv_idx.value_, new_pv_idx.value_);
    }
  }
  Assert(pv_reindexer.IsValid(GetPVCount()), "PV Reindexer is not valid.");
  return pv_reindexer;
}

// ** Explicit Instantiation
template class PartialVectorHandler<PLVType, PLVTypeEnum, NodeId>;
template class PartialVectorHandler<PLVType, PLVTypeEnum, EdgeId>;
template class PartialVectorHandler<PSVType, PSVTypeEnum, NodeId>;
template class PartialVectorHandler<PSVType, PSVTypeEnum, EdgeId>;
