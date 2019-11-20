// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A BlockModel is an abstract class to enable us to have parameter vectors that
// get subdivided and used.

#ifndef SRC_BLOCK_MODEL_HPP_
#define SRC_BLOCK_MODEL_HPP_

#include <Eigen/Dense>
#include "block_specfication.hpp"

using EigenMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using EigenVectorXdRef = Eigen::Ref<Eigen::VectorXd>;
using EigenMatrixXdRef = Eigen::Ref<EigenMatrixXd>;

class BlockModel {
 public:
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

  virtual void SetParameters(const EigenVectorXdRef parameters) = 0;

  // TODO const
  EigenVectorXdRef ExtractSegment(EigenVectorXdRef parameterization,
                                  std::string key);

  void Append(const std::string& sub_entire_key, BlockSpecification other) {
    block_specification_.Append(sub_entire_key, other);
  }

  void InsertEntireKey(BlockSpecification::Coordinates coordinates) {
    block_specification_.InsertEntireKey(coordinates);
  }

 private:
  BlockSpecification block_specification_;
  };

#endif  // SRC_BLOCK_MODEL_HPP_
