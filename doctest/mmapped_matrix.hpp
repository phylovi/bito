#pragma once

#include "../src/mmapped_matrix.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("MmappedMatrix") {
  Eigen::Index rows = 4;
  Eigen::Index cols = 5;
  using MmappedMatrixXd = MmappedMatrix<Eigen::MatrixXd>;
  {
    MmappedMatrixXd mmapped_matrix("_ignore/mmapped_matrix.data", rows, cols);
    mmapped_matrix.Get()(rows - 1, cols - 1) = 5.;
  }  // End of scope, so our mmap is destroyed and file written.
  MmappedMatrixXd mmapped_matrix("_ignore/mmapped_matrix.data", rows, cols);
  CHECK_EQ(mmapped_matrix.Get()(rows - 1, cols - 1), 5.);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
