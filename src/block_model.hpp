// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A BlockModel is an abstract class to enable us to have parameter vectors that
// get subdivided and used. To understand the structure of BlockModels, read the
// header docs for BlockSpecification first. Then have a look at the GTR model
// in substitution_model.hpp and those unit tests.
//
// The methods provided just wrap methods in BlockSpecification, so see the
// corresponding docs for a description of what they do.

#pragma once

#include "block_specification.hpp"

class BlockModel {
 public:
  BlockModel(const BlockSpecification::ParamCounts& param_counts)
      : block_specification_(param_counts) {}
  BlockModel(const BlockModel&) = delete;
  BlockModel(BlockModel&&) = delete;
  BlockModel& operator=(const BlockModel&) = delete;
  BlockModel& operator=(BlockModel&&) = delete;
  virtual ~BlockModel() = default;

  const BlockSpecification& GetBlockSpecification() const;

  EigenVectorXdRef ExtractSegment(EigenVectorXdRef param_vector, std::string key) const;
  EigenMatrixXdRef ExtractBlock(EigenMatrixXdRef param_matrix, std::string key) const;
  void Append(const std::string& sub_entire_key, BlockSpecification other);

  virtual void SetParameters(const EigenVectorXdRef param_vector) = 0;

 private:
  BlockSpecification block_specification_;
};

//  This is a virtual class so there are no unit tests. See
//  substitution_model.hpp.
