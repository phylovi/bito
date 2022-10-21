void GPEngine::SetOptimizationMethod(const GPEngine::OptimizationMethod method) {
  optimization_method_ = method;
}

void GPEngine::Optimization(const GPOperations::OptimizeBranchLength& op) {
  switch (optimization_method_) {
    case OptimizationMethod::BrentOptimization:
      return BrentOptimization(op);
    case OptimizationMethod::BrentOptimizationWithGradients:
      return BrentOptimizationWithGradients(op);
    case OptimizationMethod::GradientAscentOptimization:
      return GradientAscentOptimization(op);
    case OptimizationMethod::LogSpaceGradientAscentOptimization:
      return LogSpaceGradientAscentOptimization(op);
    case OptimizationMethod::NewtonOptimization:
      return NewtonOptimization(op);
    default:
      Failwith("GPEngine::Optimization(): Invalid OptimizationMethod given.");
  }
}

void GPEngine::SetSignificantDigitsForOptimization(int significant_digits) {
  significant_digits_for_optimization_ = significant_digits;
}

void GPEngine::BrentOptimization(const GPOperations::OptimizeBranchLength& op) {
  if (branch_length_differences_(op.gpcsp_) < branch_length_difference_threshold_) {
    return;
  }

  // Sets transition matrix to have normal branch lengths.
  auto negative_log_likelihood = [this, &op](double log_branch_length) {
    SetTransitionMatrixToHaveBranchLength(exp(log_branch_length));
    PreparePerPatternLogLikelihoodsForGPCSP(op.rootward_, op.leafward_);
    return -per_pattern_log_likelihoods_.dot(site_pattern_weights_);
  };

  double current_log_branch_length = log(branch_lengths_(op.gpcsp_));
  double current_neg_log_likelihood =
      negative_log_likelihood(current_log_branch_length);

  const auto [log_branch_length, neg_log_likelihood] =
      Optimization::BrentMinimize<false>(
          negative_log_likelihood, current_log_branch_length, min_log_branch_length_,
          max_log_branch_length_, significant_digits_for_optimization_,
          max_iter_for_optimization_, step_size_for_log_space_optimization_);

  // Numerical optimization sometimes yields new nllk > current nllk.
  // In this case, we reset the branch length to the previous value.
  if (neg_log_likelihood > current_neg_log_likelihood) {
    branch_lengths_(op.gpcsp_) = exp(current_log_branch_length);
  } else {
    branch_lengths_(op.gpcsp_) = exp(log_branch_length);
  }
  branch_length_differences_(op.gpcsp_) =
      abs(exp(current_log_branch_length) - branch_lengths_(op.gpcsp_));
}

void GPEngine::BrentOptimizationWithGradients(
    const GPOperations::OptimizeBranchLength& op) {
  if (branch_length_differences_(op.gpcsp_) < branch_length_difference_threshold_) {
    return;
  }

  auto negative_log_likelihood_and_derivative = [this, &op](double log_branch_length) {
    double branch_length = exp(log_branch_length);
    branch_lengths_(op.gpcsp_) = branch_length;
    auto [log_likelihood, log_likelihood_derivative] =
        this->LogLikelihoodAndDerivative(op);
    return std::make_pair(-log_likelihood, -branch_length * log_likelihood_derivative);
  };

  double current_log_branch_length = log(branch_lengths_(op.gpcsp_));
  double current_neg_log_likelihood =
      negative_log_likelihood_and_derivative(current_log_branch_length).first;
  const auto [log_branch_length, neg_log_likelihood] =
      Optimization::BrentMinimize<true>(
          negative_log_likelihood_and_derivative, current_log_branch_length,
          min_log_branch_length_, max_log_branch_length_,
          significant_digits_for_optimization_, max_iter_for_optimization_,
          step_size_for_log_space_optimization_);

  if (neg_log_likelihood > current_neg_log_likelihood) {
    branch_lengths_(op.gpcsp_) = exp(current_log_branch_length);
  } else {
    branch_lengths_(op.gpcsp_) = exp(log_branch_length);
  }
  branch_length_differences_(op.gpcsp_) =
      abs(exp(current_log_branch_length) - branch_lengths_(op.gpcsp_));
}

void GPEngine::GradientAscentOptimization(
    const GPOperations::OptimizeBranchLength& op) {
  auto log_likelihood_and_derivative = [this, &op](double branch_length) {
    branch_lengths_(op.gpcsp_) = branch_length;
    return this->LogLikelihoodAndDerivative(op);
  };
  const auto branch_length = Optimization::GradientAscent(
      log_likelihood_and_derivative, branch_lengths_(op.gpcsp_),
      significant_digits_for_optimization_, step_size_for_optimization_,
      min_log_branch_length_, max_iter_for_optimization_);
  branch_lengths_(op.gpcsp_) = branch_length;
}

void GPEngine::LogSpaceGradientAscentOptimization(
    const GPOperations::OptimizeBranchLength& op) {
  auto log_likelihood_and_derivative = [this, &op](double branch_length) {
    branch_lengths_(op.gpcsp_) = branch_length;
    return this->LogLikelihoodAndDerivative(op);
  };
  const auto branch_length = Optimization::LogSpaceGradientAscent(
      log_likelihood_and_derivative, branch_lengths_(op.gpcsp_),
      significant_digits_for_optimization_, step_size_for_log_space_optimization_,
      exp(min_log_branch_length_), max_iter_for_optimization_);
  branch_lengths_(op.gpcsp_) = branch_length;
}

void GPEngine::NewtonOptimization(const GPOperations::OptimizeBranchLength& op) {
  if (branch_length_differences_(op.gpcsp_) < branch_length_difference_threshold_) {
    return;
  }

  auto log_likelihood_and_first_two_derivatives = [this,
                                                   &op](double log_branch_length) {
    double x = exp(log_branch_length);
    branch_lengths_(op.gpcsp_) = x;
    auto [f_x, f_prime_x, f_double_prime_x] =
        this->LogLikelihoodAndFirstTwoDerivatives(op);
    // x = exp(y) --> f'(exp(y)) = exp(y) * f'(exp(y)) = x * f'(x)
    double f_prime_y = x * f_prime_x;
    double f_double_prime_y = f_prime_y + std::pow(x, 2) * f_double_prime_x;
    return std::make_tuple(f_x, f_prime_y, f_double_prime_y);
  };

  double current_log_branch_length = log(branch_lengths_(op.gpcsp_));
  const auto log_branch_length = Optimization::NewtonRaphsonOptimization(
      log_likelihood_and_first_two_derivatives, current_log_branch_length,
      significant_digits_for_optimization_, denominator_tolerance_for_newton_,
      min_log_branch_length_, max_log_branch_length_, max_iter_for_optimization_);

  branch_lengths_(op.gpcsp_) = exp(log_branch_length);

  branch_length_differences_(op.gpcsp_) =
      abs(exp(current_log_branch_length) - branch_lengths_(op.gpcsp_));
}
