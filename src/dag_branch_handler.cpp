// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "dag_branch_handler.hpp"

// ** Constructors

DAGBranchHandler::DAGBranchHandler(const size_t count,
                                   std::optional<OptimizationMethod> method)
    : branch_lengths_(count), differences_(count) {
  if (method.has_value()) {
    SetOptimizationMethod(method.value());
  }
  branch_lengths_.SetDefaultValue(init_default_branch_length_);
  branch_lengths_.FillWithDefault();
  differences_.SetDefaultValue(init_default_difference_);
  differences_.FillWithDefault();
}

DAGBranchHandler::DAGBranchHandler(GPDAG& dag, std::optional<OptimizationMethod> method)
    : branch_lengths_(dag), differences_(dag), dag_(&dag) {
  if (method.has_value()) {
    SetOptimizationMethod(method.value());
  }
  branch_lengths_.SetDefaultValue(init_default_branch_length_);
  branch_lengths_.FillWithDefault();
  differences_.SetDefaultValue(init_default_difference_);
  differences_.FillWithDefault();
}

// ** Comparators

int DAGBranchHandler::Compare(const DAGBranchHandler& lhs,
                              const DAGBranchHandler& rhs) {
  if (lhs.size() != rhs.size()) {
    return lhs.size() - rhs.size();
  }
  for (EdgeId edge_id{0}; edge_id < lhs.size(); edge_id++) {
    if (lhs.Get(edge_id) != rhs.Get(edge_id)) {
      return lhs.Get(edge_id) - rhs.Get(edge_id);
    }
  }
  return 0;
}

bool operator==(const DAGBranchHandler& lhs, const DAGBranchHandler& rhs) {
  return DAGBranchHandler::Compare(lhs, rhs) == 0;
}

// ** Branch Length Map

DAGBranchHandler::BranchLengthMap DAGBranchHandler::BuildBranchLengthMap(
    const GPDAG& dag) const {
  BranchLengthMap branch_length_map;
  for (EdgeId edge_id = 0; edge_id < dag.EdgeCountWithLeafSubsplits(); edge_id++) {
    auto pcsp = dag.GetDAGEdgeBitset(edge_id);
    auto branch_length = Get(edge_id);
    branch_length_map[pcsp] = branch_length;
  }
  return branch_length_map;
}

void DAGBranchHandler::ApplyBranchLengthMap(
    const DAGBranchHandler::BranchLengthMap& branch_length_map, const GPDAG& dag) {
  for (const auto& [pcsp, branch_length] : branch_length_map) {
    const EdgeId edge_id = dag.GetEdgeIdx(pcsp);
    Get(edge_id) = branch_length;
  }
}

// ** Static Functions

RootedTree DAGBranchHandler::BuildTreeWithBranchLengthsFromTopology(
    const GPDAG& dag, const DAGBranchHandler& dag_branch_handler,
    const Node::Topology& topology) {
  Tree::BranchLengthVector tree_branch_lengths(2 * topology->LeafCount() - 1, 0.0);

  topology->RootedPCSPPreorder(
      [&dag, &dag_branch_handler, &tree_branch_lengths](
          const Node* sister, const Node* focal, const Node* child_left,
          const Node* child_right) {
        Bitset parent_subsplit = Bitset::Subsplit(sister->Leaves(), focal->Leaves());
        Bitset child_subsplit =
            Bitset::Subsplit(child_left->Leaves(), child_right->Leaves());
        EdgeId focal_edge_id = dag.GetEdgeIdx(parent_subsplit, child_subsplit);
        tree_branch_lengths[focal->Id()] = dag_branch_handler(focal_edge_id);

        // If adjacent nodes go to leaves.
        if (sister->IsLeaf()) {
          Bitset subsplit = Bitset::LeafSubsplitOfNonemptyClade(sister->Leaves());
          EdgeId edge_id = dag.GetEdgeIdx(parent_subsplit, subsplit);
          tree_branch_lengths[sister->Id()] = dag_branch_handler(edge_id);
        }
        if (child_left->IsLeaf()) {
          Bitset subsplit = Bitset::LeafSubsplitOfNonemptyClade(child_left->Leaves());
          EdgeId edge_id = dag.GetEdgeIdx(child_subsplit, subsplit);
          tree_branch_lengths[child_left->Id()] = dag_branch_handler(edge_id);
        }
        if (child_right->IsLeaf()) {
          Bitset subsplit = Bitset::LeafSubsplitOfNonemptyClade(child_right->Leaves());
          EdgeId edge_id = dag.GetEdgeIdx(child_subsplit, subsplit);
          tree_branch_lengths[child_right->Id()] = dag_branch_handler(edge_id);
        }
      },
      false);

  return RootedTree(topology, std::move(tree_branch_lengths));
}

void DAGBranchHandler::CopyOverBranchLengths(const DAGBranchHandler& src,
                                             DAGBranchHandler& dest) {
  const auto& src_dag = src.GetDAG();
  const auto& dest_dag = dest.GetDAG();
  for (EdgeId dest_id = 0; dest_id < dest_dag.EdgeCountWithLeafSubsplits(); dest_id++) {
    const auto& pcsp = dest_dag.GetDAGEdgeBitset(dest_id);
    const auto src_id = src_dag.GetEdgeIdx(pcsp);
    dest.Get(dest_id) = src.Get(src_id);
  }
}

// ** Optimization

void DAGBranchHandler::OptimizeBranchLength(const EdgeId edge_id, const PVId parent_id,
                                            const PVId child_id,
                                            const bool check_branch_convergence) {
  // Branch convergence test.
  if (check_branch_convergence &&
      (differences_(edge_id) < branch_length_difference_threshold_)) {
    return;
  }

  switch (optimization_method_) {
    case OptimizationMethod::BrentOptimization:
      return BrentOptimization(edge_id, parent_id, child_id);
    case OptimizationMethod::BrentOptimizationWithGradients:
      return BrentOptimizationWithGradients(edge_id, parent_id, child_id);
    case OptimizationMethod::GradientAscentOptimization:
      return GradientAscentOptimization(edge_id, parent_id, child_id);
    case OptimizationMethod::LogSpaceGradientAscentOptimization:
      return LogSpaceGradientAscentOptimization(edge_id, parent_id, child_id);
    case OptimizationMethod::NewtonOptimization:
      return NewtonRaphsonOptimization(edge_id, parent_id, child_id);
    default:
      Failwith("DAGBranchHandler::Optimization(): Invalid OptimizationMethod given.");
  }
}

// ** Branch

void DAGBranchHandler::BrentOptimization(const EdgeId edge_id, const PVId parent_id,
                                         const PVId child_id) {
  Assert(brent_nongrad_func_ != nullptr,
         "EvalBranchFunction must be assigned before calling Brent.");
  // Evaluate branch length function.
  Optimization::NegFunc<double> negative_log_likelihood =
      [this, edge_id, parent_id, child_id](double log_branch_length) {
        return brent_nongrad_func_(edge_id, parent_id, child_id, log_branch_length);
      };
  // Capture initial.
  double current_log_branch_length = log(branch_lengths_(edge_id));
  double current_neg_log_likelihood =
      negative_log_likelihood(current_log_branch_length);
  // Optimize branch length.
  const auto [log_branch_length, neg_log_likelihood] = Optimization::BrentMinimize(
      negative_log_likelihood, current_log_branch_length, min_log_branch_length_,
      max_log_branch_length_, significant_digits_for_optimization_,
      max_iter_for_optimization_, step_size_for_log_space_optimization_);
  // Numerical optimization sometimes yields new nllk > current nllk.
  // In this case, we reset the branch length to the previous value.
  if (neg_log_likelihood > current_neg_log_likelihood) {
    branch_lengths_(edge_id) = exp(current_log_branch_length);
  } else {
    branch_lengths_(edge_id) = exp(log_branch_length);
  }
  differences_(edge_id) =
      abs(exp(current_log_branch_length) - branch_lengths_(edge_id));
}

void DAGBranchHandler::BrentOptimizationWithGradients(const EdgeId edge_id,
                                                      const PVId parent_id,
                                                      const PVId child_id) {
  Assert(brent_grad_func_ != nullptr,
         "EvalBranchFunction must be assigned before calling BrentWithGradients.");
  // Evaluate branch length function.
  Optimization::NegFuncAndDerivative<double> negative_log_likelihood_and_derivative =
      [this, edge_id, parent_id, child_id](double log_branch_length) {
        return this->brent_grad_func_(edge_id, parent_id, child_id, log_branch_length);
      };
  // Convert and capture initial branch lenth to log space.
  double current_log_branch_length = log(branch_lengths_(edge_id));
  double current_neg_log_likelihood =
      negative_log_likelihood_and_derivative(current_log_branch_length).first;
  // Optimize branch length.
  const auto [log_branch_length, neg_log_likelihood] =
      Optimization::BrentMinimizeWithGradients(
          negative_log_likelihood_and_derivative, current_log_branch_length,
          min_log_branch_length_, max_log_branch_length_,
          significant_digits_for_optimization_, max_iter_for_optimization_,
          step_size_for_log_space_optimization_);
  // Numerical optimization sometimes yields new nllk > current nllk.
  // In this case, we reset the branch length to the previous value.
  if (neg_log_likelihood > current_neg_log_likelihood) {
    branch_lengths_(edge_id) = exp(current_log_branch_length);
  } else {
    branch_lengths_(edge_id) = exp(log_branch_length);
  }
  differences_(edge_id) =
      abs(exp(current_log_branch_length) - branch_lengths_(edge_id));
}

void DAGBranchHandler::GradientAscentOptimization(const EdgeId edge_id,
                                                  const PVId parent_id,
                                                  const PVId child_id) {
  Assert(gradient_ascent_func_ != nullptr,
         "EvalBranchFunction must be assigned before calling GradientAscent.");
  // Evaluate branch length function.
  Optimization::FuncAndDerivative<double> log_likelihood_and_derivative =
      [this, edge_id, parent_id, child_id](double log_branch_length) {
        return this->gradient_ascent_func_(edge_id, parent_id, child_id,
                                           log_branch_length);
      };
  // Capture initial.
  double current_branch_length = branch_lengths_(edge_id);
  // Optimize branch length.
  const auto branch_length = Optimization::GradientAscent(
      log_likelihood_and_derivative, branch_lengths_(edge_id),
      significant_digits_for_optimization_, step_size_for_optimization_,
      min_log_branch_length_, max_iter_for_optimization_);
  // Capture result.
  branch_lengths_(edge_id) = branch_length;
  differences_(edge_id) = abs(current_branch_length - branch_lengths_(edge_id));
}

void DAGBranchHandler::LogSpaceGradientAscentOptimization(const EdgeId edge_id,
                                                          const PVId parent_id,
                                                          const PVId child_id) {
  Assert(logspace_gradient_ascent_func_ != nullptr,
         "EvalBranchFunction must be assigned before calling LogSpaceGradientAscent.");
  // Evaluate branch length function.
  Optimization::FuncAndDerivative<double> log_likelihood_and_derivative =
      [this, edge_id, parent_id, child_id](double log_branch_length) {
        return this->logspace_gradient_ascent_func_(edge_id, parent_id, child_id,
                                                    log_branch_length);
      };
  // Capture initial.
  double current_branch_length = branch_lengths_(edge_id);
  // Optimize branch length.
  const auto branch_length = Optimization::LogSpaceGradientAscent(
      log_likelihood_and_derivative, branch_lengths_(edge_id),
      significant_digits_for_optimization_, step_size_for_log_space_optimization_,
      exp(min_log_branch_length_), max_iter_for_optimization_);
  // Capture result.
  branch_lengths_(edge_id) = branch_length;
  differences_(edge_id) = abs(current_branch_length - branch_lengths_(edge_id));
}

void DAGBranchHandler::NewtonRaphsonOptimization(const EdgeId edge_id,
                                                 const PVId parent_id,
                                                 const PVId child_id) {
  Assert(newton_raphson_func_ != nullptr,
         "EvalBranchFunction must be assigned before calling NewtonRaphson.");
  // Evaluate branch length function.
  Optimization::FuncAndFirstTwoDerivatives<double>
      log_likelihood_and_first_two_derivatives =
          [this, edge_id, parent_id, child_id](double log_branch_length) {
            return this->newton_raphson_func_(edge_id, parent_id, child_id,
                                              log_branch_length);
          };
  // Capture initial.
  double current_branch_length = branch_lengths_(edge_id);
  double current_log_branch_length = log(branch_lengths_(edge_id));
  // Optimize branch length.
  const auto log_branch_length = Optimization::NewtonRaphsonOptimization(
      log_likelihood_and_first_two_derivatives, current_log_branch_length,
      significant_digits_for_optimization_, denominator_tolerance_for_newton_,
      min_log_branch_length_, max_log_branch_length_, max_iter_for_optimization_);
  // Capture result.
  branch_lengths_(edge_id) = exp(log_branch_length);
  differences_(edge_id) = abs(current_branch_length - branch_lengths_(edge_id));
}
