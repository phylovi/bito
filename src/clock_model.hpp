// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_CLOCK_MODEL_HPP_
#define SRC_CLOCK_MODEL_HPP_

#include "node.hpp"

class ClockModel {
 public:
  virtual ~ClockModel() {}

  virtual double GetRate(const Node& node) = 0;
};

class StrictClockModel : public ClockModel {
 public:
  StrictClockModel(double rate) : rate_(rate) {}

  StrictClockModel() : StrictClockModel(1.0) {}

  virtual double GetRate(const Node& node) { return rate_; }

 private:
  double rate_;
};
#endif
