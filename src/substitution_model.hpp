// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SUBSTITUTION_MODEL_HPP_
#define SRC_SUBSTITUTION_MODEL_HPP_
#include <memory>

enum struct SubstitutionModelType { JC, GTR };

class SubstitutionModel {
 public:
  virtual void EigenDecomposition() const = 0;

  static std::unique_ptr<SubstitutionModel> Create(
      SubstitutionModelType choice);
};

class GTRModel : public SubstitutionModel {
 public:
  GTRModel() = default;

  void EigenDecomposition() const override {}
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("SubstitutionModel") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SUBSTITUTION_MODEL_HPP_
