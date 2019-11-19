// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BLOCK_MODEL_HPP_
#define SRC_BLOCK_MODEL_HPP_

#include <Eigen/Dense>
#include "block_specfication.hpp"

using ParamCounts = std::unordered_map<std::string, size_t>;
using EigenRowBlock = Eigen::Ref<Eigen::VectorXd>;

class BlockModel {
 public:
  BlockModel() = default;
  BlockModel(const BlockModel&) = delete;
  BlockModel(BlockModel&&) = delete;
  BlockModel& operator=(const BlockModel&) = delete;
  BlockModel& operator=(BlockModel&&) = delete;
  virtual ~BlockModel() = default;

  virtual void SetParameters(const EigenRowBlock& parameters) = 0;

  void AddToBlockSpecification(const std::string& complete_range_name,
                               size_t& next_available_idx,
                               BlockSpecification& blocks) const;

 protected:
  virtual ParamCounts GetParamCounts() const = 0;
};

#endif  // SRC_BLOCK_MODEL_HPP_
