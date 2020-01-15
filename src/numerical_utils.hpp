//
//  numerical_utils.hpp
//  libsbn
//
//  Created by Seong-Hwan Jun on 1/13/20.
//

#ifndef SRC_NumericalUtils_HPP_
#define SRC_NumericalUtils_HPP_

#include "sugar.hpp"
#include "eigen_sugar.hpp"

namespace NumericalUtils {
  double LogAdd(double x, double y);
  double LogSum(const EigenVectorXdRef vec);
  void LogNormalize(EigenVectorXdRef vec);
  void Exponentiate(EigenVectorXdRef vec);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("numerical_utils") {
  double log_x = log(2);
  double log_y = log(3);
  double log_sum = NumericalUtils::LogAdd(log_x, log_y);
  CHECK_LT(fabs(log_sum - 1.609438), 1e-5);
  
  EigenVectorXd log_vec(10);
  double log_sum2 = DOUBLE_NEG_INF;
  for (size_t i = 0; i < log_vec.size(); i++) {
      log_vec(i) = log(i+1);
      log_sum2 = NumericalUtils::LogAdd(log_sum2, log_vec(i));
  }
  log_sum = NumericalUtils::LogSum(log_vec);
  CHECK_LT(fabs(log_sum - 4.007333), 1e-5);
  CHECK_LT(fabs(log_sum2 - 4.007333), 1e-5);
  
  // test normalization
  NumericalUtils::LogNormalize(log_vec);
  for (size_t i = 0; i < log_vec.size(); i++) {
      CHECK_LT(fabs(log_vec(i) - (log(i+1) - log_sum)), 1e-5);
  }

  // test Exponentiation
  NumericalUtils::Exponentiate(log_vec);
  double sum = 0.0;
  for (size_t i = 0; i < log_vec.size(); i++) {
    sum += log_vec(i);
  }
  CHECK_LT(fabs(sum - 1), 1e-5);
    
}
#endif  // DOCTEST_LIBRARY_INCLUDED


#endif /* SRC_NumericalUtils_HPP_ */
