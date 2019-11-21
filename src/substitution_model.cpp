// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "substitution_model.hpp"

std::unique_ptr<SubstitutionModel> SubstitutionModel::OfSpecification(
    const std::string &specification) {
  if (specification == "JC69") {
    return std::make_unique<JC69Model>();
  }  // else
  if (specification == "GTR") {
    return std::make_unique<GTRModel>();
  }  // else
  Failwith("Substitution model not known: " + specification);
}

void GTRModel::SetParameters(const EigenVectorXdRef param_vector) {
  CheckParameterVectorSize(param_vector);
  rates_ = ExtractSegment(param_vector, rates_key_);
  frequencies_ = ExtractSegment(param_vector, frequencies_key_);
  Update();
};

void GTRModel::UpdateQMatrix() {
  int rate_index = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = i + 1; j < 4; j++) {
      double rate = rates_[rate_index];
      rate_index++;
      Q_(i, j) = rate * frequencies_[j];
      Q_(j, i) = rate * frequencies_[i];
    }
  }
  // Set the diagonal entries so the rows sum to one.
  double total_substitution_rate = 0;
  for (int i = 0; i < 4; i++) {
    double row_sum = 0;
    for (int j = 0; j < 4; j++) {
      if (i != j) {
        row_sum += Q_(i, j);
      }
    }
    Q_(i, i) = -row_sum;
    total_substitution_rate += row_sum * frequencies_[i];
  }
  // Rescale matrix for unit substitution rate.
  Q_ /= total_substitution_rate;
}

void GTRModel::Update() {
  Eigen::Map<const Eigen::Array4d> tmp(&frequencies_[0]);
  EigenMatrixXd sqrt_frequencies =
      EigenMatrixXd(tmp.sqrt().matrix().asDiagonal());
  EigenMatrixXd sqrt_frequencies_inv =
      EigenMatrixXd(sqrt_frequencies.inverse());

  UpdateQMatrix();
  EigenMatrixXd S = EigenMatrixXd(sqrt_frequencies * Q_ * sqrt_frequencies_inv);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(S);

  // See p.206 of Felsenstein's book. We can get the eigendecomposition of a GTR
  // model by first getting the eigendecomposition of an associated diagonal
  // matrix and then doing this transformation.
  eigen_vectors_ = sqrt_frequencies_inv * solver.eigenvectors();
  inverse_eigen_vectors_ = solver.eigenvectors().transpose() * sqrt_frequencies;
  eigen_values_ = solver.eigenvalues();
}
