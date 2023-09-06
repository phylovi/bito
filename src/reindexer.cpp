// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "reindexer.hpp"

Reindexer Reindexer::IdentityReindexer(const size_t size) {
  Reindexer reindexer = Reindexer(size);
  std::iota(reindexer.GetData().begin(), reindexer.GetData().end(), 0);
  return reindexer;
}

bool Reindexer::IsValid(std::optional<size_t> length) const {
  const size_t reindexer_size = (length.has_value() ? length.value() : size());
  Assert(length <= size(),
         "Length of reindexer cannot be larger than the vector containing it.");
  std::vector<bool> already_used(reindexer_size, false);
  for (size_t idx = 0; idx < reindexer_size; idx++) {
    if (GetNewIndexByOldIndex(idx) >= reindexer_size ||
        already_used[GetNewIndexByOldIndex(idx)]) {
      return false;
    }
    already_used[GetNewIndexByOldIndex(idx)] = true;
  }
  return true;
}

// ** Modification Operations

void Reindexer::Resize(size_t new_size) {
  const auto old_size = GetData().size();
  GetData().resize(new_size);
  for (size_t i = old_size; i < new_size; i++) {
    GetData()[i] = i;
  }
}

// Gets the inverse of a given reindexer.
Reindexer Reindexer::InvertReindexer() const {
  Assert(IsValid(), "Reindexer must be valid in Reindexer::InvertedReindexer.");
  Reindexer inverted_reindexer(size());
  for (size_t idx = 0; idx < size(); idx++) {
    inverted_reindexer.SetReindex(GetNewIndexByOldIndex(idx), idx);
  }
  return inverted_reindexer;
}

Reindexer Reindexer::RemoveOldIndex(const size_t remove_old_idx) const {
  Assert(IsValid(), "Reindexer must be valid in Reindexer::RemoveOldIndex.");
  Reindexer result_reindexer;
  result_reindexer.reserve(size() - 1);
  const size_t remove_new_idx = GetNewIndexByOldIndex(remove_old_idx);
  for (size_t old_idx = 0; old_idx < size(); old_idx++) {
    if (old_idx == remove_old_idx) {
      continue;
    }
    const size_t new_idx = GetNewIndexByOldIndex(old_idx);
    result_reindexer.AppendNewIndex(new_idx - (new_idx > remove_new_idx));
  }
  return result_reindexer;
}

Reindexer Reindexer::RemoveNewIndex(const size_t remove_new_idx) const {
  const size_t remove_old_idx = GetOldIndexByNewIndex(remove_new_idx);
  return RemoveOldIndex(remove_old_idx);
}

Reindexer Reindexer::ComposeWith(const Reindexer &apply_reindexer) {
  Assert(IsValid(), "base_reindexer not valid in Reindexer::ComposeWith.");
  Assert(apply_reindexer.IsValid(),
         "apply_reindexer not valid in Reindexer::ComposeWith.");
  Assert(size_t(size()) <= apply_reindexer.size(),
         "The base_reindexer cannot be larger than the reindexer it is being "
         "composed with.");
  Reindexer result_reindexer(apply_reindexer.size());
  // Pad base_reindexer if it needs to grow to accept apply_reindexer.
  for (size_t idx = size(); idx < apply_reindexer.size(); idx++) {
    AppendNewIndex();
  }
  // Reindex.
  Reindexer inverted_reindexer = InvertReindexer();
  for (size_t idx = 0; idx < apply_reindexer.size(); idx++) {
    result_reindexer.SetReindex(inverted_reindexer.GetNewIndexByOldIndex(idx),
                                apply_reindexer.GetNewIndexByOldIndex(idx));
  }
  return result_reindexer;
}

void Reindexer::ReassignAndShift(const size_t old_id, const size_t new_id) {
  Assert(old_id < size() && new_id < size(),
         "The given ids must be within the bounds of the reindexer in "
         "Reindexer::ReassignAndShift.");
  Assert(IsValid(), "Reindexer must be valid in Reindexer::ReassignAndShift.");
  if (old_id == new_id) {
    return;
  }
  // Find position with value old_id.
  const size_t old_id_position = GetOldIndexByNewIndex(old_id);
  // Shift.
  if (old_id > new_id) {
    for (size_t &id : GetData()) {
      if (id < old_id && id >= new_id) {
        id++;
      }
    }
  } else {
    for (size_t &id : GetData()) {
      if (id > old_id && id <= new_id) {
        id--;
      }
    }
  }
  // Reassign old_id to new_id.
  SetReindex(old_id_position, new_id);
}
