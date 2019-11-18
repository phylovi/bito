// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "block_model.hpp"
#include "sugar.hpp"

void BlockModel::AddToBlockSpecification(BlockSpecification& blocks) const {
  size_t parameter_idx;

  for (const auto [param_name, param_length] : GetParamCounts()) {
    SafeInsert(blocks, param_name, {parameter_idx, param_length});
    parameter_idx += param_length;
  }
}
