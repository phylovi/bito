// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "pv_handler.hpp"

// ** Resize

template <class PVType, class PVTypeEnum>
void PVHandler<PVType, PVTypeEnum>::Resize(const size_t new_node_count,
                                           const size_t new_node_alloc) {
  const size_t old_plv_count = GetPLVCount();
  node_count_ = new_node_count;
  node_alloc_ = new_node_alloc;
  // Allocate mmapped data block.
  mmapped_master_plvs_.Resize(GetAllocatedPLVCount() * pattern_count_);
  // Subdivide mmapped data in individual PVs.
  plvs_ = mmapped_master_plvs_.Subdivide(GetAllocatedPLVCount());
  // Initialize new work space.
  Assert((plvs_.back().rows() == MmappedNucleotidePLV::base_count_) &&
             (plvs_.back().cols() == static_cast<Eigen::Index>(pattern_count_)) &&
             (size_t(plvs_.size()) == GetAllocatedPLVCount()),
         "Didn't get the right shape of PVs out of Subdivide.");
  for (size_t i = old_plv_count; i < GetPaddedPLVCount(); i++) {
    plvs_.at(i).setZero();
  }
}

template <class PVType, class PVTypeEnum>
void PVHandler<PVType, PVTypeEnum>::Reindex(const Reindexer plv_reindexer) {
  Reindexer::ReindexInPlace(plvs_, plv_reindexer, GetPLVCount(), GetPLV(GetPLVCount()),
                            GetPLV(GetPLVCount() + 1));
}

template <class PVType, class PVTypeEnum>
Reindexer PVHandler<PVType, PVTypeEnum>::BuildPLVReindexer(
    const Reindexer &node_reindexer, const size_t old_node_count,
    const size_t new_node_count) {
  node_count_ = new_node_count;
  Reindexer plv_reindexer(new_node_count * plv_count_per_node_);
  size_t new_plvs_idx = old_node_count * plv_count_per_node_;
  for (size_t i = 0; i < new_node_count; i++) {
    const size_t old_node_idx = i;
    const size_t new_node_idx = node_reindexer.GetNewIndexByOldIndex(old_node_idx);
    for (const auto plv_type : typename PVTypeEnum::Iterator()) {
      // Either get input plv_index from old plvs, or get new plv_index (new data is
      // irrelevant, so just get next available index).
      size_t old_plv_idx;
      if (old_node_idx < old_node_count) {
        old_plv_idx = GetPLVIndex(plv_type, old_node_idx, old_node_count);
      } else {
        old_plv_idx = new_plvs_idx;
        new_plvs_idx++;
      }
      const size_t new_plv_idx = GetPLVIndex(plv_type, new_node_idx, new_node_count);
      plv_reindexer.SetReindex(old_plv_idx, new_plv_idx);
    }
  }
  Assert(plv_reindexer.IsValid(GetPLVCount()), "PV Reindexer is not valid.");
  return plv_reindexer;
}

// explicit instatiation
template class PVHandler<PartialVectorType::PLVType, PartialVectorType::PLVTypeEnum>;
template class PVHandler<PartialVectorType::PSVType, PartialVectorType::PSVTypeEnum>;
