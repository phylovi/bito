#pragma once

#include "../src/substitution_model.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <algorithm>
TEST_CASE("SubstitutionModel") {
  auto CheckEigenvalueEquality = [](EigenVectorXd eval1, EigenVectorXd eval2) {
    std::sort(eval1.begin(), eval1.end());
    std::sort(eval2.begin(), eval2.end());
    CheckVectorXdEquality(eval1, eval2, 0.0001);
  };
  auto gtr_model = std::make_unique<GTRModel>();
  auto hky_model = std::make_unique<HKYModel>();
  auto jc_model = std::make_unique<JC69Model>();
  // Test 1: First we test using the "built in" default values.
  CheckEigenvalueEquality(jc_model->GetEigenvalues(), gtr_model->GetEigenvalues());
  CheckEigenvalueEquality(jc_model->GetEigenvalues(), hky_model->GetEigenvalues());
  EigenVectorXd param_vector(10);
  // Test 2: Now try out ParameterSegmentMapOf.
  gtr_model = std::make_unique<GTRModel>();
  // First zero out our param_vector.
  param_vector.setZero();
  // We can use ParameterSegmentMapOf to get two "views" into our parameter
  // vector.
  auto parameter_map =
      gtr_model->GetBlockSpecification().ParameterSegmentMapOf(param_vector);
  auto frequencies = parameter_map.at(SubstitutionModel::frequencies_key_);
  auto rates = parameter_map.at(SubstitutionModel::rates_key_);
  // When we modify the contents of these views, that changes param_vector.
  frequencies.setConstant(0.25);
  rates.setConstant(1.0 / 6.0);
  // We can then set param_vector and go forward as before.
  gtr_model->SetParameters(param_vector);
  CheckEigenvalueEquality(jc_model->GetEigenvalues(), gtr_model->GetEigenvalues());
  // Test 3: Compare to eigenvalues from R.
  frequencies << 0.479367, 0.172572, 0.140933, 0.207128;
  rates << 0.060602, 0.402732, 0.028230, 0.047910, 0.407249, 0.053277;
  gtr_model->SetParameters(param_vector);
  EigenVectorXd eigen_values_r(4);
  eigen_values_r << -2.567992e+00, -1.760838e+00, -4.214918e-01, 1.665335e-16;
  CheckEigenvalueEquality(eigen_values_r, gtr_model->GetEigenvalues());
  // Test HKY against GTR
  EigenVectorXd hky_param_vector(5);
  hky_param_vector.setZero();
  auto hky_parameter_map =
      hky_model->GetBlockSpecification().ParameterSegmentMapOf(hky_param_vector);
  auto hky_frequencies = hky_parameter_map.at(SubstitutionModel::frequencies_key_);
  auto hky_kappa = hky_parameter_map.at(SubstitutionModel::rates_key_);
  hky_frequencies << 0.1, 0.2, 0.3, 0.4;
  hky_kappa.setConstant(3.0);
  hky_model->SetParameters(hky_param_vector);
  frequencies << 0.1, 0.2, 0.3, 0.4;
  rates << 0.1, 0.3, 0.1, 0.1, 0.3, 0.1;
  gtr_model->SetParameters(param_vector);
  CheckEigenvalueEquality(gtr_model->GetEigenvalues(), hky_model->GetEigenvalues());
  CHECK(gtr_model->GetQMatrix().isApprox(hky_model->GetQMatrix()));
}
#endif  // DOCTEST_LIBRARY_INCLUDED
