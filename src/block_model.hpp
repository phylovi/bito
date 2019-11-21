// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A BlockModel is an abstract class to enable us to have parameter vectors that
// get subdivided and used. To understand the structure of BlockModels, read the
// header docs for BlockSpecification first. Then have a look at the GTR model
// in substitution_model.hpp and those unit tests.

#ifndef SRC_BLOCK_MODEL_HPP_
#define SRC_BLOCK_MODEL_HPP_

#include <Eigen/Dense>
#include "block_specification.hpp"

using EigenMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using EigenVectorXdRef = Eigen::Ref<Eigen::VectorXd>;
using EigenMatrixXdRef = Eigen::Ref<EigenMatrixXd>;

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

  const BlockSpecification& GetBlockSpecification() const {
    return block_specification_;
  }

  void CheckParametersSize(const EigenVectorXdRef parameters) const {
    Assert(parameters.size() == block_specification_.ParameterCount(),
           "Parameters are the wrong dimension!");
  }
  void CheckParameterMatrixSize(const EigenMatrixXdRef parameter_matrix) const {
    Assert(parameter_matrix.cols() == block_specification_.ParameterCount(),
           "Parameters are the wrong dimension!");
  }

  virtual void SetParameters(const EigenVectorXdRef parameters) = 0;

  EigenVectorXdRef ExtractSegment(EigenVectorXdRef parameterization,
                                  std::string key) const;
  EigenMatrixXdRef ExtractBlock(EigenMatrixXdRef parameter_matrix,
                                std::string key) const;

  ParameterSegmentMap ParameterSegmentMapOf(EigenVectorXdRef parameters) const;
  ParameterBlockMap ParameterBlockMapOf(
      EigenMatrixXdRef parameter_matrix) const;

  void Append(const std::string& sub_entire_key, BlockSpecification other) {
    block_specification_.Append(sub_entire_key, other);
  }

  // TODO comment out.
  void InsertEntireKey(BlockSpecification::Coordinates coordinates) {
    block_specification_.InsertEntireKey(coordinates);
  }

 private:
  BlockSpecification block_specification_;
  };

#endif  // SRC_BLOCK_MODEL_HPP_
