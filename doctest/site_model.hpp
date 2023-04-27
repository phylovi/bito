#pragma once

#include "../src/site_model.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <algorithm>
TEST_CASE("SiteModel") {
  // Test 1: First we test using the "built in" default values.
  auto weibull_model = std::make_unique<WeibullSiteModel>(4, 1.0);
  const EigenVectorXd rates = weibull_model->GetCategoryRates();
  EigenVectorXd rates_r(4);
  rates_r << 0.1457844, 0.5131316, 1.0708310, 2.2702530;
  CheckVectorXdEquality(rates, rates_r, 0.0001);

  // Test 2: Now set param_vector using SetParameters.
  weibull_model = std::make_unique<WeibullSiteModel>(4, 1.0);
  EigenVectorXd param_vector(1);
  param_vector << 0.1;
  weibull_model->SetParameters(param_vector);
  rates_r << 4.766392e-12, 1.391131e-06, 2.179165e-03, 3.997819e+00;
  const EigenVectorXd rates2 = weibull_model->GetCategoryRates();
  CheckVectorXdEquality(rates2, rates_r, 0.0001);

  // Test 3: Check proportions.
  const EigenVectorXd proportions = weibull_model->GetCategoryProportions();
  CheckVectorXdEquality(0.25, proportions, 0.0001);

  // Test 4: Check sum rates[i]*proportions[i]==1.
  CHECK_LT(fabs(rates.dot(proportions) - 1.), 0.0001);
  CHECK_LT(fabs(rates2.dot(proportions) - 1.), 0.0001);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
