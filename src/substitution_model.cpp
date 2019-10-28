// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "substitution_model.hpp"

void GTRModel::UpdateEigenDecomposition() {
  Eigen::Map<const Eigen::Array4d> tmp(&frequencies_[0]);
  EigenMatrix4d sqrt_frequencies = tmp.sqrt().matrix().asDiagonal();
  EigenMatrix4d sqrt_frequencies_inv = sqrt_frequencies.inverse();
  EigenMatrix4d Q;

  int rate_index = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = i + 1; j < 4; j++) {
      double rate = rates_[rate_index];
      rate_index++;
      Q(i, j) = rate * frequencies_[j];
      Q(j, i) = rate * frequencies_[i];
    }
  }
  // Set the diagonal entries so the rows sum to one.
  double total_substitution_rate = 0;
  for (int i = 0; i < 4; i++) {
    double row_sum = 0;
    for (int j = 0; j < 4; j++) {
      if (i != j) {
        row_sum += Q(i, j);
      }
    }
    Q(i, i) = -row_sum;
    total_substitution_rate += row_sum * frequencies_[i];
  }
  // Rescale matrix for unit substitution rate.
  Q /= total_substitution_rate;

  EigenMatrix4d S = EigenMatrix4d(sqrt_frequencies * Q * sqrt_frequencies_inv);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(S);

  const EigenMatrix4d& eigen_vectors =
      sqrt_frequencies_inv * solver.eigenvectors();
  const EigenMatrix4d& inverse_eigen_vectors =
      solver.eigenvectors().transpose() * sqrt_frequencies;
  const EigenVector4d& eigen_values = solver.eigenvalues();

  // TODO Did you try using the eigen STL iterators?
  // https://eigen.tuxfamily.org/dox-devel/group__TutorialSTL.html
  std::copy(&eigen_vectors.data()[0], &eigen_vectors.data()[0] + 16,
            evec_.begin());
  std::copy(&inverse_eigen_vectors.data()[0],
            &inverse_eigen_vectors.data()[0] + 16, ivec_.begin());
  std::copy(&eigen_values.data()[0], &eigen_values.data()[0] + 4,
            eval_.begin());

  need_update_ = false;
}
