// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SUBSTITUTION_MODEL_HPP_
#define SRC_SUBSTITUTION_MODEL_HPP_
#include <memory>

struct SubstitutionModelDetails {
  size_t state_count;
};

class SubstitutionModel {
 public:
  virtual void EigenDecomposition() const = 0;
  virtual SubstitutionModelDetails Details() const = 0;
};
using SubstitutionModelPtr = std::unique_ptr<SubstitutionModel>;

class GTRModel : public SubstitutionModel {
 public:
  GTRModel() : details_({4}) {}

  SubstitutionModelDetails Details() const override { return details_; }

  void EigenDecomposition() const override {}

 private:
  SubstitutionModelDetails details_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("SubstitutionModel") {
  auto substitution_model = std::make_unique<GTRModel>();

  CHECK_EQ(substitution_model->Details().state_count, 4);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SUBSTITUTION_MODEL_HPP_
