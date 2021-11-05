// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <memory>
#include <string>

#include "block_model.hpp"
#include "node.hpp"

class ClockModel : public BlockModel {
 public:
  ClockModel(const BlockSpecification::ParamCounts& param_counts)
      : BlockModel(param_counts) {}
  virtual ~ClockModel() = default;

  virtual double GetRate(size_t node_id) = 0;

  static std::unique_ptr<ClockModel> OfSpecification(const std::string& specification);
};

class NoClockModel : public ClockModel {
 public:
  explicit NoClockModel() : ClockModel({}) {}

  double GetRate(size_t node_id) override { return 1.; }

  void SetParameters(const EigenVectorXdRef parameters) override{};
};

class StrictClockModel : public ClockModel {
 public:
  explicit StrictClockModel(double rate) : ClockModel({{rate_key_, 1}}), rate_(rate) {}

  StrictClockModel() : StrictClockModel(1.0) {}

  double GetRate(size_t node_id) override { return rate_; }

  void SetParameters(const EigenVectorXdRef parameters) override;

  inline const static std::string rate_key_ = "clock rate";

 private:
  double rate_;
};
