// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SUBSTITUTION_MODEL_HPP_
#define SRC_SUBSTITUTION_MODEL_HPP_
#include <memory>

enum struct SubstitutionModelType { JC, GTR };

struct SubstitutionModelDetails {
  SubstitutionModelType type;
};

class SubstitutionModel {
 public:
  SubstitutionModelDetails Details() const { return details_; }

  virtual void EigenDecomposition() const = 0;

 protected:
  explicit SubstitutionModel(SubstitutionModelDetails details)
      : details_(details) {}
  SubstitutionModelDetails details_;
};
using SubstitutionModelPtr = std::unique_ptr<SubstitutionModel>;

class JCModel : public SubstitutionModel {
 public:
  JCModel() : SubstitutionModel({SubstitutionModelType::JC}) {}

  void EigenDecomposition() const override {}
};

class GTRModel : public SubstitutionModel {
 public:
  GTRModel() : SubstitutionModel({SubstitutionModelType::GTR}) {}

  void EigenDecomposition() const override {}
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("SubstitutionModel") {
  auto substitution_model = std::make_unique<GTRModel>();

  CHECK_EQ(substitution_model->Details().type, SubstitutionModelType::GTR);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SUBSTITUTION_MODEL_HPP_
