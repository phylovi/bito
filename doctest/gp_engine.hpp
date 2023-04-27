#pragma once

#include "../src/gp_engine.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("GPEngine") {
  EigenVectorXd empty_vector;
  SitePattern hello_site_pattern = SitePattern::HelloSitePattern();
  GPEngine engine(hello_site_pattern, 5, 5, "_ignore/mmapped_plv.data",
                  GPEngine::default_rescaling_threshold_, empty_vector, empty_vector,
                  empty_vector, false);
  engine.SetTransitionMatrixToHaveBranchLength(0.75);
  // Computed directly:
  // https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_%28Jukes_and_Cantor_1969%29
  CHECK(fabs(0.52590958087 - engine.GetTransitionMatrix()(0, 0)) < 1e-10);
  CHECK(fabs(0.1580301397 - engine.GetTransitionMatrix()(0, 1)) < 1e-10);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
