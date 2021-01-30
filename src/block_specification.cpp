// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "block_specification.hpp"

void BlockSpecification::Insert(const std::string& key, Coordinates value) {
  SafeInsert(map_, key, value);
}

void BlockSpecification::Insert(const char* key, Coordinates value) {
  Insert(std::string(key), value);
}

BlockSpecification::BlockSpecification(ParamCounts param_counts) {
  size_t next_available_idx = 0;
  for (const auto [block_name, block_size] : param_counts) {
    Insert(block_name, {next_available_idx, block_size});
    next_available_idx += block_size;
  }
  InsertEntireKey({0, next_available_idx});
}

BlockSpecification::Coordinates BlockSpecification::Find(const std::string& key) const {
  auto search = map_.find(key);
  if (search == map_.end()) {
    Failwith("Can't find '" + key + "' in block specification!");
  }
  return search->second;
}

void BlockSpecification::InsertEntireKey(Coordinates coordinates) {
  EraseEntireKey();
  Insert(entire_key_, coordinates);
}

void BlockSpecification::Append(const std::string& sub_entire_key,
                                BlockSpecification other) {
  const auto our_parameter_count = ParameterCount();
  for (const auto [block_name, coordinate] : other.GetMap()) {
    auto [start_idx, block_size] = coordinate;
    if (block_name == entire_key_) {
      Assert(start_idx == 0, "Start index of entire block isn't zero.");
      Insert(sub_entire_key, {our_parameter_count, block_size});
    } else {
      auto search = map_.find(block_name);
      if (search != map_.end()) {
        Failwith("Key overlap between BlockSpecifications: " + block_name);
      }  // else
      Insert(block_name, {our_parameter_count + start_idx, block_size});
    }
  }
  InsertEntireKey({0, our_parameter_count + other.ParameterCount()});
}

void BlockSpecification::CheckParameterVectorSize(
    const EigenVectorXdRef param_vector) const {
  Assert(param_vector.size() == ParameterCount(),
         "Parameters are the wrong dimension!");
}
void BlockSpecification::CheckParameterMatrixSize(
    const EigenMatrixXdRef param_matrix) const {
  Assert(param_matrix.cols() == ParameterCount(),
         "Parameters are the wrong dimension!");
}

EigenVectorXdRef BlockSpecification::ExtractSegment(EigenVectorXdRef param_vector,
                                                    std::string key) const {
  auto [start_idx, parameter_count] = Find(key);
  if (start_idx + parameter_count > param_vector.size()) {
    Failwith("Model parameter '" + key +
             "' request too long for a param_vector of length " +
             std::to_string(param_vector.size()) + ".");
  }  // else
  return param_vector.segment(start_idx, parameter_count);
}

EigenMatrixXdRef BlockSpecification::ExtractBlock(EigenMatrixXdRef param_matrix,
                                                  std::string key) const {
  auto [start_idx, parameter_count] = Find(key);
  if (start_idx + parameter_count > param_matrix.cols()) {
    Failwith("Model parameter '" + key +
             "' request too long for a param_matrix of width " +
             std::to_string(param_matrix.cols()) + ".");
  }  // else
  return param_matrix.block(0, start_idx, param_matrix.rows(), parameter_count);
}

// In the interest of having a single code path, this implementation does a
// lookups in ExtractSegment even though it doesn't really need to. My intent is
// not to have it be called frequently-- just once when we're setting things up.
// So if this changes let's do something about it.
BlockSpecification::ParameterSegmentMap BlockSpecification::ParameterSegmentMapOf(
    EigenVectorXdRef param_vector) const {
  CheckParameterVectorSize(param_vector);
  ParameterSegmentMap parameter_segment_map;
  for (const auto [key, _] : GetMap()) {
    SafeInsert(parameter_segment_map, key, ExtractSegment(param_vector, key));
  }
  return parameter_segment_map;
}

BlockSpecification::ParameterBlockMap BlockSpecification::ParameterBlockMapOf(
    EigenMatrixXdRef param_matrix) const {
  ParameterBlockMap parameter_block_map;
  CheckParameterMatrixSize(param_matrix);
  for (const auto [key, _] : GetMap()) {
    SafeInsert(parameter_block_map, key, ExtractBlock(param_matrix, key));
  }
  return parameter_block_map;
}
