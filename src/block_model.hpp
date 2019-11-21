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
  // These are handy structures that turn the block specification into a map
  // from the block specification keys to segments (i.e. sub-vectors) and blocks
  // (i.e. sub-matrices). We can then read and write to these values, which will
  // be reflected in the original parameter vector/matrix.
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

  void CheckParameterVectorSize(const EigenVectorXdRef param_vector) const;
  void CheckParameterMatrixSize(const EigenMatrixXdRef param_matrix) const;

  // These methods allow us to pull out segments (i.e. sub-vectors) from vectors
  // and blocks from matrices depending on the coordinates of the BlockModel.
  // They are very useful for writing the SetParameters method.
  EigenVectorXdRef ExtractSegment(EigenVectorXdRef param_vector,
                                  std::string key) const;
  EigenMatrixXdRef ExtractBlock(EigenMatrixXdRef param_matrix,
                                std::string key) const;

  // These are explained in the definition of ParameterSegmentMap and
  // ParameterBlockMap.
  ParameterSegmentMap ParameterSegmentMapOf(
      EigenVectorXdRef param_vector) const;
  ParameterBlockMap ParameterBlockMapOf(
      EigenMatrixXdRef param_matrix) const;

  void Append(const std::string& sub_entire_key, BlockSpecification other) {
    block_specification_.Append(sub_entire_key, other);
  }

  virtual void SetParameters(const EigenVectorXdRef param_vector) = 0;

 private:
  BlockSpecification block_specification_;
  };

  //  This is a virtual class so there are no unit tests. See
  //  substitution_model.hpp.
#endif  // SRC_BLOCK_MODEL_HPP_
