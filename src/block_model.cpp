// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "block_model.hpp"
#include "sugar.hpp"

void BlockModel::AddToBlockSpecification(size_t& next_available_idx,
                                         BlockSpecification& blocks) const {
  for (const auto [param_name, param_length] : GetParamCounts()) {
    SafeInsert(blocks, param_name, {next_available_idx, param_length});
    next_available_idx += param_length;
  }
}
