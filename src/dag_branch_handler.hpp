// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This class handles data on the DAG in the form of a DAGData object and DAG branch
// lengths, including storage and resizing of data vectors and branch length
// optimization.

#pragma once

#include "gp_dag.hpp"
#include "dag_data.hpp"
#include "optimization.hpp"
#include "substitution_model.hpp"
#include "tree.hpp"

class DAGBranchHandler {
 public:
  DAGBranchHandler(const size_t count,
                   std::optional<OptimizationMethod> method = std::nullopt);
  DAGBranchHandler(GPDAG& dag, std::optional<OptimizationMethod> method = std::nullopt);

  // ** Comparators

  static int Compare(const DAGBranchHandler& lhs, const DAGBranchHandler& rhs);
  friend bool operator==(const DAGBranchHandler& lhs, const DAGBranchHandler& rhs);

  // ** Counts

  size_t GetCount() { return GetBranchLengths().GetCount(); }
  size_t GetSpareCount() { return GetBranchLengths().GetSpareCount(); }
  size_t GetAllocCount() { return GetBranchLengths().GetAllocCount(); }

  void SetCount(const size_t count) {
    GetBranchLengths().SetCount(count);
    GetBranchDifferences().SetCount(count);
  }
  void SetSpareCount(const size_t count) {
    GetBranchLengths().SetSpareCount(count);
    GetBranchDifferences().SetSpareCount(count);
  }
  void SetAllocCount(const size_t count) {
    GetBranchLengths().SetAllocCount(count);
    GetBranchDifferences().SetAllocCount(count);
  }

  // Number of optimization iterations called.
  size_t GetOptimizationCount() const { return optimization_count_; }
  // Reset optimization iterations to zero, and reset branch length differences to zero.
  void ResetOptimizationCount() {
    optimization_count_ = 0;
    differences_.FillWithDefault();
  }
  // Add one optimization iteration to count.
  void IncrementOptimizationCount() { optimization_count_++; }
  // Checks if called current optimization iteration is the first.
  bool IsFirstOptimization() const { return optimization_count_ == 0; }

  // ** Access

  // Get the size of the data vector.
  size_t size() const { return branch_lengths_.size(); }

  // Get underlying DAGData objects for Branch Lengths.
  DAGEdgeDoubleData& GetBranchLengths() { return branch_lengths_; }
  const DAGEdgeDoubleData& GetBranchLengths() const { return branch_lengths_; }
  void SetBranchLengths(EigenVectorXd& branch_lengths) {
    GetBranchLengths().SetDAGData(branch_lengths);
  }
  // Get underlying DAGData objects for Branch Differences.
  DAGEdgeDoubleData& GetBranchDifferences() { return differences_; }
  const DAGEdgeDoubleData& GetBranchDifferences() const { return differences_; }
  void SetBranchDifferences(EigenVectorXd& differences) {
    GetBranchDifferences().SetDAGData(differences);
  }

  // Get underlying EigenVector for Branch Lengths.
  EigenVectorXd& GetBranchLengthData() { return branch_lengths_.GetData(); }
  const EigenVectorXd& GetBranchLengthData() const { return branch_lengths_.GetData(); }
  // Get underlying EigenVector for Branch Differences.
  EigenVectorXd& GetBranchDifferenceData() { return differences_.GetData(); }
  const EigenVectorXd& GetBranchDifferenceData() const {
    return differences_.GetData();
  }

  // Access elements in data vector.
  double& Get(const EdgeId edge_id) { return branch_lengths_(edge_id); }
  const double& Get(const EdgeId edge_id) const { return branch_lengths_(edge_id); }
  double& operator()(const EdgeId edge_id) { return Get(edge_id); }
  const double& operator()(const EdgeId edge_id) const { return Get(edge_id); }
  double& GetBranchDifference(const EdgeId edge_id) { return differences_(edge_id); }
  const double& GetBranchDifference(const EdgeId edge_id) const {
    return differences_(edge_id);
  }

  // References.
  const bool HasDAG() const { return dag_ != nullptr; }
  const void SetDAG(GPDAG& dag) { dag_ = &dag; }
  const GPDAG& GetDAG() const {
    Assert(HasDAG(), "Reference DAG cannot be accessed before it has been set.");
    return *dag_;
  }
  const bool HasGraftDAG() const { return graft_dag_ != nullptr; }
  const void SetGraftDAG(GraftDAG& graft_dag) { graft_dag_ = &graft_dag; }
  const GraftDAG& GetGraftDAG() const {
    Assert(HasGraftDAG(),
           "Reference GraftDAG cannot be accessed before it has been set.");
    return *graft_dag_;
  }

  // Method used for branch length optimization.
  void SetOptimizationMethod(const OptimizationMethod method) {
    optimization_method_ = method;
  }
  OptimizationMethod GetOptimizationMethod() const { return optimization_method_; }
  // Set number of significant digits of precision used in branch length optimization.
  void SetSignificantDigitsForOptimization(int significant_digits) {
    significant_digits_for_optimization_ = significant_digits;
  }
  double GetSignificantDigitsForOptimization() const {
    return significant_digits_for_optimization_;
  }
  // Set initial values for branch length optimization.
  void SetDefaultBranchLength(double default_branch_length) {
    branch_lengths_.SetDefaultValue(default_branch_length);
  }
  double GetDefaultBranchLength() const { return branch_lengths_.GetDefaultValue(); }

  // ** Resize

  // Resize using given count.
  void Resize(std::optional<const size_t> edge_count = std::nullopt,
              std::optional<const size_t> spare_count = std::nullopt,
              std::optional<const size_t> explicit_alloc = std::nullopt,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    branch_lengths_.Resize(edge_count, spare_count, explicit_alloc, reindexer);
    differences_.Resize(edge_count, spare_count, explicit_alloc, reindexer);
  }
  // Resize using reference DAG.
  void Resize(std::optional<const size_t> explicit_alloc = std::nullopt,
              std::optional<const Reindexer> reindexer = std::nullopt) {
    branch_lengths_.Resize(GetDAG(), explicit_alloc, reindexer);
    differences_.Resize(GetDAG(), explicit_alloc, reindexer);
  }

  // ** Evaluation Functions
  // Each optimization method needs an evaluation function should take in an edge, its
  // parent and child nodes, and a branch length, and return a set of specified values
  // according to the optimization method, such as the log likelihood and its
  // derivatives.
  // - Takes the following args: (1) focal EdgeId, (2) rootward PVId (parent node's
  // RFocal PV), (3) leafward PVId (child node's P PV), (4) log branch length of focal
  // edge.
  using LogLikelihoodAndDerivativeFunc =
      std::function<DoublePair(EdgeId, PVId, PVId, double)>;
  using LogLikelihoodAndFirstTwoDerivativesFunc =
      std::function<Tuple3<double>(EdgeId, PVId, PVId, double)>;
  using NegLogLikelihoodFunc = std::function<double(EdgeId, PVId, PVId, double)>;
  using NegLogLikelihoodAndDerivativeFunc =
      std::function<DoublePair(EdgeId, PVId, PVId, double)>;

  // Set helper function for Nongradient Brent. Function takes a branch length and
  // returns the negative log likelihood.
  void SetBrentFunc(NegLogLikelihoodFunc brent_nongrad_func) {
    brent_nongrad_func_ = brent_nongrad_func;
  }
  NegLogLikelihoodFunc GetBrentFunc() { return brent_nongrad_func_; }
  // Set helper function for Gradient Brent. Function takes a branch length and returns
  // the negative log likelihood and the first derivative.
  void SetBrentWithGradientFunc(NegLogLikelihoodAndDerivativeFunc brent_grad_func) {
    brent_grad_func_ = brent_grad_func;
  }
  NegLogLikelihoodAndDerivativeFunc GetBrentWithGradientFunc() {
    return brent_grad_func_;
  }
  // Set helper function for Gradient Ascent. Function should return the log likelihood
  // and the first derivative.
  void SetGradientAscentFunc(LogLikelihoodAndDerivativeFunc gradient_ascent_func) {
    gradient_ascent_func_ = gradient_ascent_func;
  }
  LogLikelihoodAndDerivativeFunc GetGradientAscentFunc() {
    return gradient_ascent_func_;
  }
  // Set helper function for Gradient Ascent. Function takes a log branch length and
  // returns the log likelihood and the first derivative.
  void SetLogSpaceGradientAscentFunc(
      LogLikelihoodAndDerivativeFunc logspace_gradient_ascent_func) {
    logspace_gradient_ascent_func_ = logspace_gradient_ascent_func;
  }
  LogLikelihoodAndDerivativeFunc GetLogSpaceGradientAscentFunc() {
    return logspace_gradient_ascent_func_;
  }
  // Set helper function for Newton-Raphson. Function takes a log branch length and
  // returns the log likelihood and the first two derivatives.
  void SetNewtonRaphsonFunc(
      LogLikelihoodAndFirstTwoDerivativesFunc newton_raphson_func) {
    newton_raphson_func_ = newton_raphson_func;
  }
  LogLikelihoodAndFirstTwoDerivativesFunc GetNewtonRaphsonFunc() {
    return newton_raphson_func_;
  }

  // ** Branch Length Optimization

  // Performs optimization on given branch. Parent_id expects the parent node's RFocal
  // PV. Child_id expects the child node's P PV.
  void OptimizeBranchLength(const EdgeId edge_id, const PVId parent_rfocal_pvid,
                            const PVId child_p_pvid,
                            const bool check_branch_convergence = true);

  // ** Branch Length Map

  // Build a map from PCSP bitset to branch length.
  using BranchLengthMap = std::map<Bitset, double>;
  BranchLengthMap BuildBranchLengthMap(const GPDAG& dag) const;

  // Copy over branch lengths from map by corresponding PCSPs. Only applies branch
  // lengths for PCSPs which occur in DAG, all others left unchanged.  All edges in map
  // must exist in dag.
  void ApplyBranchLengthMap(const BranchLengthMap& branch_length_map, const GPDAG& dag);

  // ** Static Functions

  // Builds a tree using given branch lengths on a given topology. For each edge, builds
  // out PCSP bitset to find corresponding DAG EdgeId, to find branch length.
  static RootedTree BuildTreeWithBranchLengthsFromTopology(
      const GPDAG& dag, const DAGBranchHandler& dag_branch_handler,
      const Node::Topology& topology);

  // Copies branch lengths from one handler to another. Base DAG of dest handler must be
  // a subgraph of src handler.
  static void CopyOverBranchLengths(const DAGBranchHandler& src,
                                    DAGBranchHandler& dest);

 protected:
  // ** Branch Length Optimization Helpers

  void BrentOptimization(const EdgeId edge_id, const PVId parent_id,
                         const PVId child_id);
  void BrentOptimizationWithGradients(const EdgeId edge_id, const PVId parent_id,
                                      const PVId child_id);
  void GradientAscentOptimization(const EdgeId edge_id, const PVId parent_id,
                                  const PVId child_id);
  void LogSpaceGradientAscentOptimization(const EdgeId edge_id, const PVId parent_id,
                                          const PVId child_id);
  void NewtonRaphsonOptimization(const EdgeId edge_id, const PVId parent_id,
                                 const PVId child_id);

  // Branch Lengths.
  DAGEdgeDoubleData branch_lengths_;
  // Tracks differences between branch length over different iterations of
  // optimization to test for edge-wise convergence.
  DAGEdgeDoubleData differences_;

  // Unowned DAG reference.
  GPDAG* dag_ = nullptr;
  // Unowned GraftDAG reference.
  GraftDAG* graft_dag_ = nullptr;

  // Optimization pass counter.  Tracks number of iterations of optimization.
  size_t optimization_count_ = 0;
  // Method used in branch length optimization.
  OptimizationMethod optimization_method_ = OptimizationMethod::BrentOptimization;

  // Default branch length set at initialization. Current branch lengths are stored
  // in branch_lengths_.
  static constexpr double init_default_branch_length_ = 0.1;
  // Default difference set at initialization. Current branch lengths are stored
  // in differences_.
  static constexpr double init_default_difference_ = 0.0;
  // Absolute lower bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double min_log_branch_length_ = -13.9;
  // Absolute upper bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double max_log_branch_length_ = 1.1;
  // Precision used for checking convergence of branch length optimization.
  // In the non-Brent optimization methods, significant digits will be used to
  // determine convergence through relative tolerance, i.e. measuring difference
  // from previous branch length values until the absolute difference is below
  // 10^(-significant_digits_for_optimization_).
  // Brent optimization does not define convergence through relative tolerance,
  // rather convergence based on tightness of the brackets that it adapts during
  // optimization. This variable thus represents the "number of bits precision to
  // which the minima should be found". When testing on sample datasets, we found that
  // setting the value to 10 was a good compromise between speed and precision for
  // Brent. See more on Brent optimization here:
  // https://www.boost.org/doc/libs/1_79_0/libs/math/doc/html/math_toolkit/brent_minima.html
  int significant_digits_for_optimization_ = 10;
  double relative_tolerance_for_optimization_ = 1e-4;
  double denominator_tolerance_for_newton_ = 1e-10;
  double step_size_for_optimization_ = 5e-4;
  double step_size_for_log_space_optimization_ = 1.0005;
  // Number of iterations allowed for branch length optimization.
  size_t max_iter_for_optimization_ = 1000;
  double branch_length_difference_threshold_ = 1e-15;

  // Evaluation Functions
  NegLogLikelihoodFunc brent_nongrad_func_ = nullptr;
  NegLogLikelihoodAndDerivativeFunc brent_grad_func_ = nullptr;
  LogLikelihoodAndDerivativeFunc gradient_ascent_func_ = nullptr;
  LogLikelihoodAndDerivativeFunc logspace_gradient_ascent_func_ = nullptr;
  LogLikelihoodAndFirstTwoDerivativesFunc newton_raphson_func_ = nullptr;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

#endif  // DOCTEST_LIBRARY_INCLUDED
