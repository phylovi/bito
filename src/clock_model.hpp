// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_CLOCK_MODEL_HPP_
#define SRC_CLOCK_MODEL_HPP_

#include "block_model.hpp"
#include "node.hpp"

class ClockModel : public BlockModel {
 public:
  ClockModel() = default;
  virtual ~ClockModel() = default;

  virtual double GetRate(const Node& node) = 0;

  static std::unique_ptr<ClockModel> OfSpecification(
      const std::string& specification);
  void SetParameters(const EigenVectorRef parameters){};
};

class StrictClockModel : public ClockModel {
 public:
  explicit StrictClockModel(double rate) : rate_(rate) {}

  StrictClockModel() : StrictClockModel(1.0) {}

  double GetRate(const Node& node) override { return rate_; }

  ParamCounts GetParamCounts() const override { return {{"clock rate", 1}}; }

 private:
  double rate_;
};
#endif
