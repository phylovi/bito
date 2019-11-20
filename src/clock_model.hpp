// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_CLOCK_MODEL_HPP_
#define SRC_CLOCK_MODEL_HPP_

#include "block_model.hpp"
#include "node.hpp"

class ClockModel : public BlockModel {
 public:
  ClockModel(const BlockSpecification::ParamCounts& param_counts)
      : BlockModel(param_counts) {}
  virtual ~ClockModel() = default;

  virtual double GetRate(const Node& node) = 0;

  static std::unique_ptr<ClockModel> OfSpecification(
      const std::string& specification);
};

class StrictClockModel : public ClockModel {
 public:
  explicit StrictClockModel(double rate)
      : ClockModel({{rate_key_, 1}}), rate_(rate) {}

  StrictClockModel() : StrictClockModel(1.0) {}

  double GetRate(const Node& node) override { return rate_; }

  void SetParameters(const EigenVectorXdRef parameters) override {
    Assert(parameters.size() == 1,
           "StrictClockModel parameters are the wrong length!");
    rate_ = parameters[0];
  };

  inline const static std::string rate_key_ = "clock rate";

 private:
  double rate_;
};
#endif
