// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "reindexer.hpp"

Reindexer Reindexer::IdentityReindexer(const size_t size) {
  Reindexer reindexer = Reindexer(size);
  reindexer.IdentityReindexer();
  return reindexer;
}

// ** Access

size_t Reindexer::GetInputIndexByOutputIndex(
    const size_t new_index, std::optional<Reindexer> inverted_reindexer) const {
  if (inverted_reindexer.has_value()) {
    return inverted_reindexer.value().GetOutputIndexByInputIndex(new_index);
  }
  return size_t(std::find(GetData().begin(), GetData().end(), new_index) -
                GetData().begin());
}

// ** Modify

void Reindexer::ReassignOutputIndexAndShift(const size_t old_output_idx,
                                            const size_t new_output_idx) {
  Assert(old_output_idx < size() && new_output_idx < size(),
         "The given ids must be within the bounds of the reindexer in "
         "Reindexer::ReassignOutputIndexAndShift.");
  Assert(IsValid(),
         "Reindexer must be valid in Reindexer::ReassignOutputIndexAndShift.");
  if (old_output_idx == new_output_idx) {
    return;
  }
  // Find position with value old_output_idx.
  const size_t old_input_idx = GetInputIndexByOutputIndex(old_output_idx);
  // Shift.
  if (old_output_idx > new_output_idx) {
    for (size_t &id : GetData()) {
      if (id < old_output_idx && id >= new_output_idx) {
        id++;
      }
    }
  } else {
    for (size_t &id : GetData()) {
      if (id > old_output_idx && id <= new_output_idx) {
        id--;
      }
    }
  }
  // Reassign old_output_idx to new_output_idx.
  SetReindex(old_input_idx, new_output_idx);
}

void Reindexer::AppendNextIndex(const size_t append_count) {
  for (size_t i = 0; i < append_count; i++) {
    data_.push_back(size());
  }
}

void Reindexer::AppendOutputIndex(const size_t new_index) {
  data_.push_back(new_index);
}

void Reindexer::RemoveInputIndex(const size_t input_idx_to_remove) {
  const size_t output_idx_to_remove = GetOutputIndexByInputIndex(input_idx_to_remove);
  for (size_t idx = 0; idx < size(); idx++) {
    // Skip index we are removing.
    if (idx == input_idx_to_remove) {
      continue;
    }
    // Close gap for old and new indices if we are past the removed index.
    const size_t input_idx = (idx - (input_idx > input_idx_to_remove));
    const size_t output_idx = GetOutputIndexByInputIndex(input_idx);
    data_[input_idx] = (output_idx - (output_idx > output_idx_to_remove));
  }
  data_.pop_back();
}

void Reindexer::RemoveOutputIndex(const size_t output_idx_to_remove) {
  const size_t input_idx_to_remove = GetInputIndexByOutputIndex(output_idx_to_remove);
  std::cout << "input_idx->output_idx_to_remove: " << input_idx_to_remove << " "
            << output_idx_to_remove << std::endl;
  return RemoveInputIndex(input_idx_to_remove);
}

void Reindexer::RemoveInputIndex(SizeVector &input_idx_to_remove) {
  SizeVector output_idx_to_remove(input_idx_to_remove);
  for (size_t i = 0; i < input_idx_to_remove.size(); i++) {
    output_idx_to_remove[i] = GetOutputIndexByInputIndex(input_idx_to_remove[i]);
  }
  std::sort(input_idx_to_remove.begin(), input_idx_to_remove.end());
  input_idx_to_remove.push_back(data_.size());
  std::sort(output_idx_to_remove.begin(), output_idx_to_remove.end());
  output_idx_to_remove.push_back(data_.size());
  // Close gaps create by deletions.
  for (size_t i = 0; i < input_idx_to_remove.size() - 1; i++) {
    const size_t padding = i + 1;
    for (size_t j = input_idx_to_remove[i] + 1; j < input_idx_to_remove[i + 1]; j++) {
      size_t input_idx = j;
      size_t output_idx = data_[input_idx];
      for (size_t k = 0; k < output_idx_to_remove.size() - 1; k++) {
        if (output_idx < output_idx_to_remove[k]) {
          break;
        }
        output_idx--;
      }
      data_[input_idx - padding] = output_idx;
    }
  }
  for (size_t i = 0; i < input_idx_to_remove.size() - 1; i++) {
    data_.pop_back();
  }
}

void Reindexer::RemoveOutputIndex(SizeVector &output_idx_to_remove,
                                  std::optional<Reindexer> inverted_reindexer) {
  // Convert all new indices to old indices, then call other routine.
  SizeVector input_idx_to_remove;
  if (!inverted_reindexer.has_value()) {
    inverted_reindexer = InvertReindexer();
  }
  for (size_t i = 0; i < output_idx_to_remove.size(); i++) {
    input_idx_to_remove[i] =
        GetInputIndexByOutputIndex(output_idx_to_remove[i], inverted_reindexer);
  }
  return RemoveInputIndex(input_idx_to_remove);
}

void Reindexer::InsertInputIndex(const size_t input_idx_to_add) {
  SizeVector input_idx_to_add_vec({input_idx_to_add});
  InsertInputIndex(input_idx_to_add_vec);
}

void Reindexer::InsertOutputIndex(const size_t output_idx_to_add) {
  SizeVector output_idx_to_add_vec({output_idx_to_add});
  InsertOutputIndex(output_idx_to_add_vec);
}

void Reindexer::InsertOutputIndex(SizeVector &sorted_output_idxs_to_add) {
  if (sorted_output_idxs_to_add.size() == 0) {
    return;
  }
  const size_t old_size = data_.size();
  // Padding vector gives spans for different levels of padding.
  SizeVector padding_vector(sorted_output_idxs_to_add);
  padding_vector.push_back(data_.size());
  std::sort(padding_vector.begin(), padding_vector.end());
  // Allocate space for new indices.
  AppendNextIndex(sorted_output_idxs_to_add.size());
  // Pad out space for new indices.
  for (size_t i = padding_vector.size() - 2; i >= 1; i--) {
    const size_t padding = i + 1;
    for (size_t j = padding_vector[i + 1] - 1; j >= padding_vector[i]; j--) {
      data_[j + padding] = data_[j];
    }
  }
  // Insert new values.
  for (size_t i = 0; i < sorted_output_idxs_to_add.size(); i++) {
    data_[sorted_output_idxs_to_add[i]] = old_size + i;
  }
}

void Reindexer::InsertInputIndex(SizeVector &sorted_input_idxs_to_add) {
  if (sorted_input_idxs_to_add.size() == 0) {
    return;
  }
  const size_t old_size = data_.size();
  // pad to account for previous inserted indices
  for (size_t i = 0; i < sorted_input_idxs_to_add.size(); i++) {
    sorted_input_idxs_to_add[i] = sorted_input_idxs_to_add[i] + i;
  }
  // Append new indexes to end.
  data_.insert(data_.end(), sorted_input_idxs_to_add.begin(),
               sorted_input_idxs_to_add.end());
  // Shift indices to accomodate new indices.
  for (size_t i = 0; i < old_size; i++) {
    const size_t data_idx = GetOutputIndexByInputIndex(i);
    for (size_t j = 0; j < sorted_input_idxs_to_add.size(); j++) {
      const size_t inserted_idx = sorted_input_idxs_to_add[j];
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
    inverted_reindexer.SetReindex(GetOutputIndexByInputIndex(idx), idx);
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
    result_reindexer.SetReindex(inverted_reindexer.GetOutputIndexByInputIndex(idx),
                                apply_reindexer.GetOutputIndexByInputIndex(idx));
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
    if (GetOutputIndexByInputIndex(idx) >= reindexer_size ||
        already_used[GetOutputIndexByInputIndex(idx)]) {
      return false;
    }
    already_used[GetOutputIndexByInputIndex(idx)] = true;
  }
  return true;
}
