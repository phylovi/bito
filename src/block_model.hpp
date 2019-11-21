// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A BlockModel is an abstract class to enable us to have parameter vectors that
// get subdivided and used. To understand the structure of BlockModels, read the
// header docs for BlockSpecification first. Then have a look at the GTR model
// in substitution_model.hpp and those unit tests.

#ifndef SRC_BLOCK_MODEL_HPP_
#define SRC_BLOCK_MODEL_HPP_

#include "block_specification.hpp"
#include "eigen_sugar.hpp"

class BlockModel {
 public:
  using ParameterSegmentMap = std::unordered_map<std::string, EigenVectorXdRef>;
  using ParameterBlockMap = std::unordered_map<std::string, EigenMatrixXdRef>;

  BlockModel(const BlockSpecification::ParamCounts& param_counts)
      : block_specification_(param_counts) {}
  BlockModel(const BlockModel&) = delete;
  BlockModel(BlockModel&&) = delete;
  BlockModel& operator=(const BlockModel&) = delete;
  BlockModel& operator=(BlockModel&&) = delete;
  virtual ~BlockModel() = default;

  const BlockSpecification& GetBlockSpecification() const;

  virtual void SetParameters(const EigenVectorXdRef parameters) = 0;

  void CheckParametersSize(const EigenVectorXdRef parameters) const;
  void CheckParameterMatrixSize(const EigenMatrixXdRef param_matrix) const;

  EigenVectorXdRef ExtractSegment(EigenVectorXdRef parameterization,
                                  std::string key) const;
  EigenMatrixXdRef ExtractBlock(EigenMatrixXdRef param_matrix,
                                std::string key) const;

  ParameterSegmentMap ParameterSegmentMapOf(EigenVectorXdRef parameters) const;
  ParameterBlockMap ParameterBlockMapOf(
      EigenMatrixXdRef param_matrix) const;

  void Append(const std::string& sub_entire_key, BlockSpecification other) {
    block_specification_.Append(sub_entire_key, other);
  }

 private:
  BlockSpecification block_specification_;
  };

  //  This is a virtual class so no unit tests. See substitution_model.hpp.
#endif  // SRC_BLOCK_MODEL_HPP_
