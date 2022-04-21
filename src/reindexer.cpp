// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "reindexer.hpp"

Reindexer Reindexer::IdentityReindexer(const size_t size) {
  Reindexer reindexer = Reindexer(size);
  reindexer.IdentityReindexer();
  return reindexer;
}

// ** Access

size_t Reindexer::GetOldIndexByNewIndex(
    const size_t new_index, std::optional<Reindexer> inverted_reindexer) const {
  if (inverted_reindexer.has_value()) {
    return inverted_reindexer.value().GetNewIndexByOldIndex(new_index);
  }
  return size_t(std::find(GetData().begin(), GetData().end(), new_index) -
                GetData().begin());
}

// ** Modify

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

void Reindexer::AppendNextIndex(const size_t append_count) {
  for (size_t i = 0; i < append_count; i++) {
    data_.push_back(size());
  }
}

void Reindexer::AppendNewIndex(const size_t new_index) { data_.push_back(new_index); }

void Reindexer::RemoveOldIndex(const size_t old_idx_to_remove) {
  const size_t new_idx_to_remove = GetNewIndexByOldIndex(old_idx_to_remove);
  for (size_t idx = 0; idx < size(); idx++) {
    // Skip index we are removing.
    if (idx == old_idx_to_remove) {
      continue;
    }
    // Close gap for old and new indices if we are past the removed index.
    const size_t old_idx = (idx - (old_idx > old_idx_to_remove));
    const size_t new_idx = GetNewIndexByOldIndex(old_idx);
    data_[old_idx] = (new_idx - (new_idx > new_idx_to_remove));
  }
  data_.pop_back();
}

void Reindexer::RemoveNewIndex(const size_t new_idx_to_remove) {
  const size_t old_idx_to_remove = GetOldIndexByNewIndex(new_idx_to_remove);
  return RemoveOldIndex(old_idx_to_remove);
}

void Reindexer::RemoveOldIndex(SizeVector &old_idx_to_remove) {}

void Reindexer::RemoveNewIndex(SizeVector &new_idx_to_remove,
                               std::optional<Reindexer> inverted_reindexer) {
  // Convert all new indices to old indices, then call other routine.
  SizeVector old_idx_to_remove;
  if (!inverted_reindexer.has_value()) {
    inverted_reindexer = InvertReindexer();
  }
  for (size_t i = 0; i < new_idx_to_remove.size(); i++) {
    old_idx_to_remove[i] =
        GetOldIndexByNewIndex(new_idx_to_remove[i], inverted_reindexer);
  }
  return RemoveOldIndex(old_idx_to_remove);
}

void Reindexer::InsertOldIndex(const size_t old_idx_to_add) {
  SizeVector old_idx_to_add_vec({old_idx_to_add});
  InsertOldIndex(old_idx_to_add_vec);
}

void Reindexer::InsertNewIndex(const size_t new_idx_to_add) {
  SizeVector new_idx_to_add_vec({new_idx_to_add});
  InsertNewIndex(new_idx_to_add_vec);
}

void Reindexer::InsertNewIndex(SizeVector &sorted_new_idxs_to_add) {
  if (sorted_new_idxs_to_add.size() == 0) {
    return;
  }
  const size_t old_size = data_.size();
  // Padding vector gives spans for different levels of padding.
  SizeVector padding_vector(sorted_new_idxs_to_add);
  padding_vector.push_back(data_.size());
  std::sort(padding_vector.begin(), padding_vector.end());
  // Allocate space for new indices.
  AppendNextIndex(sorted_new_idxs_to_add.size());
  // Pad out space for new indices.
  for (size_t i = padding_vector.size() - 2; i >= 1; i--) {
    size_t padding = i + 1;
    for (size_t j = padding_vector[i + 1] - 1; j >= padding_vector[i]; j--) {
      data_[j + padding] = data_[j];
    }
  }
  // Insert new values.
  for (size_t i = 0; i < sorted_new_idxs_to_add.size(); i++) {
    data_[sorted_new_idxs_to_add[i]] = old_size + i;
  }
}

void Reindexer::InsertOldIndex(SizeVector &sorted_old_idxs_to_add) {
  if (sorted_old_idxs_to_add.size() == 0) {
    return;
  }
  const size_t old_size = data_.size();
  // pad to account for previous inserted indices
  for (size_t i = 0; i < sorted_old_idxs_to_add.size(); i++) {
    sorted_old_idxs_to_add[i] = sorted_old_idxs_to_add[i] + i;
  }
  // Append new indexes to end.
  data_.insert(data_.end(), sorted_old_idxs_to_add.begin(),
               sorted_old_idxs_to_add.end());
  // Shift indices to accomodate new indices.
  for (size_t i = 0; i < old_size; i++) {
    const size_t data_idx = GetNewIndexByOldIndex(i);
    for (size_t j = 0; j < sorted_old_idxs_to_add.size(); j++) {
      const size_t inserted_idx = sorted_old_idxs_to_add[j];
      if (data_idx < inserted_idx) {
        break;
      }
      data_[i]++;
    }
  }
}

// ** Transform

void Reindexer::IdentityReindexer() {
  std::iota(GetData().begin(), GetData().end(), 0);
}

Reindexer Reindexer::InvertReindexer() const {
  Assert(IsValid(), "Reindexer must be valid in Reindexer::InvertReindexer.");
  Reindexer inverted_reindexer(size());
  for (size_t idx = 0; idx < size(); idx++) {
    inverted_reindexer.SetReindex(GetNewIndexByOldIndex(idx), idx);
  }
  return inverted_reindexer;
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
    AppendNextIndex();
  }
  // Reindex.
  Reindexer inverted_reindexer = InvertReindexer();
  for (size_t idx = 0; idx < apply_reindexer.size(); idx++) {
    result_reindexer.SetReindex(inverted_reindexer.GetNewIndexByOldIndex(idx),
                                apply_reindexer.GetNewIndexByOldIndex(idx));
  }
  return result_reindexer;
}

// ** Miscellaneous

std::ostream &operator<<(std::ostream &os, const Reindexer &reindexer) {
  os << reindexer.GetData();
  return os;
};

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
