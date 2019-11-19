// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "block_model.hpp"
#include "sugar.hpp"

void BlockModel::AddToBlockSpecification(const std::string& complete_range_name,
                                         size_t& next_available_idx,
                                         BlockSpecification& blocks) const {
  size_t original_next_available_idx = next_available_idx;
  size_t total_param_length = 0;
  for (const auto [param_name, param_length] : GetParamCounts()) {
    blocks.Insert(param_name, {next_available_idx, param_length});
    next_available_idx += param_length;
    total_param_length += param_length;
  }
  blocks.Insert(complete_range_name,
                {original_next_available_idx, total_param_length});
}
