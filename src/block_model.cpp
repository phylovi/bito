// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "block_model.hpp"
#include "sugar.hpp"

const BlockSpecification& BlockModel::GetBlockSpecification() const {
  return block_specification_;
}

void BlockModel::CheckParameterVectorSize(
    const EigenVectorXdRef param_vector) const {
  Assert(param_vector.size() == block_specification_.ParameterCount(),
         "Parameters are the wrong dimension!");
}
void BlockModel::CheckParameterMatrixSize(
    const EigenMatrixXdRef param_matrix) const {
  Assert(param_matrix.cols() == block_specification_.ParameterCount(),
         "Parameters are the wrong dimension!");
}

EigenVectorXdRef BlockModel::ExtractSegment(EigenVectorXdRef param_vector,
                                            std::string key) const {
  auto [start_idx, parameter_count] = block_specification_.Find(key);
  if (start_idx + parameter_count > param_vector.size()) {
    Failwith("Model parameter " + key +
             " request too long for a param_vector of length " +
             std::to_string(param_vector.size()) + ".");
  }  // else
  return param_vector.segment(start_idx, parameter_count);
}

EigenMatrixXdRef BlockModel::ExtractBlock(EigenMatrixXdRef param_matrix,
                                          std::string key) const {
  auto [start_idx, parameter_count] = block_specification_.Find(key);
  if (start_idx + parameter_count > param_matrix.cols()) {
    Failwith("Model parameter " + key +
             " request too long for a param_matrix of width " +
             std::to_string(param_matrix.cols()) + ".");
  }  // else
  return param_matrix.block(0, start_idx, param_matrix.rows(),
                                parameter_count);
}

// In the interest of having a single code path, this implementation does twice
// as many lookups as it needs by calling ExtractSegment. My intent is not to
// have it be called frequently-- just once when we're setting things up. So if
// this changes let's do something about it.
BlockModel::ParameterSegmentMap BlockModel::ParameterSegmentMapOf(
    EigenVectorXdRef param_vector) const {
  CheckParameterVectorSize(param_vector);
  ParameterSegmentMap parameter_segment_map;
  for (const auto [key, _] : GetBlockSpecification().GetMap()) {
    SafeInsert(parameter_segment_map, key,
               ExtractSegment(param_vector, key));
  }
  return parameter_segment_map;
}

BlockModel::ParameterBlockMap BlockModel::ParameterBlockMapOf(
    EigenMatrixXdRef param_matrix) const {
  ParameterBlockMap parameter_block_map;
  CheckParameterMatrixSize(param_matrix);
  for (const auto [key, _] : GetBlockSpecification().GetMap()) {
    SafeInsert(parameter_block_map, key, ExtractBlock(param_matrix, key));
  }
  return parameter_block_map;
}
