// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BLOCK_MODEL_HPP_
#define SRC_BLOCK_MODEL_HPP_

#include <Eigen/Dense>
#include <unordered_map>
#include "sugar.hpp"

using BlockCoordinates = std::tuple<size_t, size_t>;
using BlockSpecification = std::unordered_map<std::string, BlockCoordinates>;
using ParamCounts = std::unordered_map<std::string, size_t>;

class BlockModel {
 public:
  // TODO delete default constructor?
  BlockModel() = default;
  explicit BlockModel(size_t parameter_start)
      : parameter_count_(0), parameter_start_(parameter_start){};
  BlockModel(const BlockModel&) = delete;
  BlockModel(BlockModel&&) = delete;
  BlockModel& operator=(const BlockModel&) = delete;
  BlockModel& operator=(BlockModel&&) = delete;
  virtual ~BlockModel() = default;

  virtual void SetParameters(const Eigen::VectorXd& parameters) = 0;

  void AddToBlockSpecification(size_t& next_available_idx,
                               BlockSpecification& blocks) const;

 protected:
  size_t parameter_count_;
  size_t parameter_start_;

  virtual ParamCounts GetParamCounts() const = 0;
};

#endif  // SRC_BLOCK_MODEL_HPP_
