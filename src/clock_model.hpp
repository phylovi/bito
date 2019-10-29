// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_CLOCK_MODEL_HPP_
#define SRC_CLOCK_MODEL_HPP_

#include "node.hpp"

class ClockModel {
 public:
  // TODO Will we need a custom destructor here? I'm suspecting that we'll just
  // have vectors and such as member variables, which will just go away ("rule
  // of zero").
  virtual ~ClockModel() = default;

  virtual double GetRate(const Node& node) = 0;
};

class StrictClockModel : public ClockModel {
 public:
  explicit StrictClockModel(double rate) : rate_(rate) {}

  StrictClockModel() : StrictClockModel(1.0) {}

  double GetRate(const Node& node) override { return rate_; }

 private:
  double rate_;
};
#endif
