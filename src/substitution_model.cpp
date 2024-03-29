// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "substitution_model.hpp"

std::unique_ptr<SubstitutionModel> SubstitutionModel::OfSpecification(
    const std::string &specification) {
  if (specification == "JC69") {
    return std::make_unique<JC69Model>();
  }  // else
  if (specification == "HKY") {
    return std::make_unique<HKYModel>();
  }  // else
  if (specification == "GTR") {
    return std::make_unique<GTRModel>();
  }  // else
  Failwith("Substitution model not known: " + specification);
}

void JC69Model::UpdateEigendecomposition() {
  eigenvectors_ << 1.0, 2.0, 0.0, 0.5, 1.0, -2.0, 0.5, 0.0, 1.0, 2.0, 0.0, -0.5, 1.0,
      -2.0, -0.5, 0.0;
  inverse_eigenvectors_ << 0.25, 0.25, 0.25, 0.25, 0.125, -0.125, 0.125, -0.125, 0.0,
      1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0;
  eigenvalues_ << 0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333;
}

void JC69Model::UpdateQMatrix() {
  Q_ << -1.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, -1.0, 1.0 / 3.0, 1.0 / 3.0,
      1.0 / 3.0, 1.0 / 3.0, -1.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, -1.0;
}

void HKYModel::SetParameters(const EigenVectorXdRef param_vector) {
  GetBlockSpecification().CheckParameterVectorSize(param_vector);
  rates_ = ExtractSegment(param_vector, rates_key_);
  frequencies_ = ExtractSegment(param_vector, frequencies_key_);
  if (fabs(frequencies_.sum() - 1.) >= 0.001) {
    std::ostringstream oss;
    std::copy(frequencies_.begin(), frequencies_.end() - 1,
              std::ostream_iterator<double>(oss, ","));
    oss << *frequencies_.end();
    Failwith("HKY frequencies do not sum to 1 +/- 0.001! frequency vector: (" +
             oss.str() + ")");
  }
  Update();
};

void HKYModel::UpdateQMatrix() {
  EigenVectorXd rates;
  rates.resize(6);
  rates << 1.0, rates_[0], 1.0, 1.0, rates_[0], 1.0;
  int rate_index = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = i + 1; j < 4; j++) {
      double rate = rates[rate_index];
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

// Analytical eigendecomposition.
// See Hasegawa, Kishino, and Yano, 1985 for details.
void HKYModel::UpdateEigendecomposition() {
  double kappa = rates_[0];

  double pi_a = frequencies_[0];
  double pi_c = frequencies_[1];
  double pi_g = frequencies_[2];
  double pi_t = frequencies_[3];

  double pi_r = pi_a + pi_g;
  double pi_y = pi_c + pi_t;

  double beta = -1.0 / (2.0 * (pi_r * pi_y + kappa * (pi_a * pi_g + pi_c * pi_t)));
  eigenvalues_[0] = 0;
  eigenvalues_[1] = beta;
  eigenvalues_[2] = beta * (1 + pi_y * (kappa - 1));
  eigenvalues_[3] = beta * (1 + pi_r * (kappa - 1));

  inverse_eigenvectors_.setZero();
  eigenvectors_.setZero();

  inverse_eigenvectors_.row(0) << pi_a, pi_c, pi_g, pi_t;

  inverse_eigenvectors_.row(1) << pi_a * pi_y, -pi_c * pi_r, pi_g * pi_y, -pi_t * pi_r;

  inverse_eigenvectors_(2, 1) = 1;
  inverse_eigenvectors_(2, 3) = -1;

  inverse_eigenvectors_(3, 0) = 1;
  inverse_eigenvectors_(3, 2) = -1;

  eigenvectors_.col(0).setOnes();

  eigenvectors_.col(1) << 1. / pi_r, -1. / pi_y, 1. / pi_r, -1. / pi_y;

  eigenvectors_(1, 2) = pi_t / pi_y;
  eigenvectors_(3, 2) = -pi_c / pi_y;

  eigenvectors_(0, 3) = pi_g / pi_r;
  eigenvectors_(2, 3) = -pi_a / pi_r;
}

void GTRModel::SetParameters(const EigenVectorXdRef param_vector) {
  GetBlockSpecification().CheckParameterVectorSize(param_vector);
  rates_ = ExtractSegment(param_vector, rates_key_);
  frequencies_ = ExtractSegment(param_vector, frequencies_key_);
  if (fabs(frequencies_.sum() - 1.) >= 0.001) {
    std::ostringstream oss;
    std::copy(frequencies_.begin(), frequencies_.end() - 1,
              std::ostream_iterator<double>(oss, ","));
    oss << *frequencies_.end();
    Failwith("GTR frequencies do not sum to 1 +/- 0.001! frequency vector: (" +
             oss.str() + ")");
  }
  if (fabs(rates_.sum() - 1.) >= 0.001) {
    std::ostringstream oss;
    std::copy(rates_.begin(), rates_.end() - 1,
              std::ostream_iterator<double>(oss, ","));
    oss << *rates_.end();
    Failwith("GTR rates do not sum to 1 +/- 0.001! rate vector: (" + oss.str() + ")");
  }
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

void DNAModel::UpdateEigendecomposition() {
  Eigen::Map<const Eigen::Array4d> tmp(&frequencies_[0]);
  EigenMatrixXd sqrt_frequencies = EigenMatrixXd(tmp.sqrt().matrix().asDiagonal());
  EigenMatrixXd sqrt_frequencies_inv = EigenMatrixXd(sqrt_frequencies.inverse());

  EigenMatrixXd S = EigenMatrixXd(sqrt_frequencies * Q_ * sqrt_frequencies_inv);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(S);

  // See p.206 of Felsenstein's book. We can get the eigendecomposition of a GTR
  // model by first getting the eigendecomposition of an associated diagonal
  // matrix and then doing this transformation.
  eigenvectors_ = sqrt_frequencies_inv * solver.eigenvectors();
  inverse_eigenvectors_ = solver.eigenvectors().transpose() * sqrt_frequencies;
  eigenvalues_ = solver.eigenvalues();
}

void DNAModel::Update() {
  UpdateQMatrix();
  UpdateEigendecomposition();
}
