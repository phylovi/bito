// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "block_model.hpp"

#include "sugar.hpp"

const BlockSpecification& BlockModel::GetBlockSpecification() const {
  return block_specification_;
}

EigenVectorXdRef BlockModel::ExtractSegment(EigenVectorXdRef param_vector,
                                            std::string key) const {
  return block_specification_.ExtractSegment(param_vector, key);
}

EigenMatrixXdRef BlockModel::ExtractBlock(EigenMatrixXdRef param_matrix,
                                          std::string key) const {
  return block_specification_.ExtractBlock(param_matrix, key);
}

void BlockModel::Append(const std::string& sub_entire_key, BlockSpecification other) {
  block_specification_.Append(sub_entire_key, other);
}
