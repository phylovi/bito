// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// TP Evaluation Engine is an interface for the TPEngine, that facilitates different
// methods for scoring Top Trees. Each edge in the DAG corresponds to its "top tree",
// which is the best scoring tree contained in the DAG has the given edge, which is
// stored in TPEngine's choice maps. There are various methods for evaluating tree
// scores.  Currently tree likelihoods and tree parsimony are supported.
//
// Note: Alongside scoring, handles DAG data particular to given methods. For example,
// the likelhood engine maintains DAG branch lengths.

#pragma once

#include "sugar.hpp"
#include "gp_dag.hpp"
#include "graft_dag.hpp"
#include "pv_handler.hpp"
#include "tp_choice_map.hpp"
#include "nni_operation.hpp"
#include "sankoff_handler.hpp"
#include "dag_branch_handler.hpp"
#include "optimization.hpp"
#include "substitution_model.hpp"
#include "stopwatch.hpp"

class TPEngine;
using BitsetEdgeIdMap = std::unordered_map<Bitset, EdgeId>;

// TPEngine helper for evaluating Top Trees.
class TPEvalEngine {
 public:
  TPEvalEngine(TPEngine &tp_engine);

  // ** Maintenance

  // Initialize Computation Engine.
  virtual void Initialize() = 0;
  // Update the Computation Engine after adding Node Pairs to the DAG.
  virtual void UpdateEngineAfterDAGAddNodePair(const NNIOperation &post_nni,
                                               const NNIOperation &pre_nni,
                                               std::optional<size_t> new_tree_id) = 0;
  // Update the Computation Engine after modifying DAG. Nodes and Edges in reindexer
  // that exceed the prev count gives the location of the new node and edge ids.
  virtual void UpdateEngineAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer) = 0;
  // Computes scores. Call after Initialize or any Update steps, and before
  // GetTopTreeScores.
  virtual void ComputeScores() = 0;

  // ** Scoring

  // Get the Top Tree from the DAG with the given edge.
  virtual double GetTopTreeScoreWithEdge(const EdgeId edge_id) const;
  // Get the Top Tree from the DAG containing the proposed NNI.
  virtual double GetTopTreeScoreWithProposedNNI(
      const NNIOperation &post_nni, const NNIOperation &pre_nni,
      const size_t spare_offset = 0,
      std::optional<BitsetEdgeIdMap> best_edge_map = std::nullopt) = 0;

  // ** Resize

  // Resize Engine for modified DAG.
  virtual void GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                                std::optional<Reindexer> edge_reindexer);
  // Grow engine to handle computing NNIs for all adjacent NNIs.
  // Option to grow engine for computing via reference or via copy. If computing via
  // reference, option whether to use unique temporaries (for testing and computing in
  // parallel).
  virtual void GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                         const bool via_reference = true,
                                         const bool use_unique_temps = true);

  // Resize GPEngine to accomodate DAG with given number of nodes and edges.  Option
  // to remap data according to DAG reindexers.  Option to give explicit number of
  // nodes or edges to allocate memory for (this is the only way memory allocation
  // will be decreased).
  virtual void GrowNodeData(
      const size_t node_count,
      std::optional<const Reindexer> node_reindexer = std::nullopt,
      std::optional<const size_t> explicit_alloc = std::nullopt,
      const bool on_init = false);
  virtual void GrowEdgeData(
      const size_t edge_count,
      std::optional<const Reindexer> edge_reindexer = std::nullopt,
      std::optional<const size_t> explicit_alloc = std::nullopt,
      const bool on_init = false);
  // Grow space for storing temporary computation.
  virtual void GrowSpareNodeData(const size_t new_node_spare_count);
  virtual void GrowSpareEdgeData(const size_t new_edge_spare_count);

  // Copy all edge data from its pre_edge_id to post_edge_id.
  virtual void CopyEdgeData(const EdgeId src_edge_id, const EdgeId dest_edge_id);

  // ** Access

  // Get reference NNIEngine.
  TPEngine &GetTPEngine() { return *tp_engine_; }
  const TPEngine &GetTPEngine() const { return *tp_engine_; }
  // Get reference DAG.
  const GPDAG &GetDAG() const { return *dag_; }
  // Get reference GraftDAG.
  const GraftDAG &GetGraftDAG() const { return *graft_dag_; }
  // Get reference SitePattern.
  const SitePattern &GetSitePattern() const { return *site_pattern_; }
  // Get top tree score vector by edge.
  EigenVectorXd &GetTopTreeScores() { return top_tree_per_edge_; }
  const EigenVectorXd &GetTopTreeScores() const { return top_tree_per_edge_; }

 protected:
  // Un-owned reference to NNIEngine.
  TPEngine *tp_engine_ = nullptr;
  // Un-owned reference to DAG.
  const GPDAG *dag_ = nullptr;
  // Un-owned reference to GraftDAG.
  const GraftDAG *graft_dag_ = nullptr;
  // Observed leave states.
  SitePattern *site_pattern_ = nullptr;

  // Top-scoring tree per edge.
  EigenVectorXd top_tree_per_edge_;
};

// TPEngine helper for evaluating Top Trees using Likelihood.
class TPEvalEngineViaLikelihood : public TPEvalEngine {
 public:
  TPEvalEngineViaLikelihood(TPEngine &tp_engine, const std::string &mmap_path);

  // ** Maintenance

  // Initialize Computation Engine and Populate PVs.
  void Initialize() override;
  // Update the Computation Engine after adding Node Pairs to the DAG.
  void UpdateEngineAfterDAGAddNodePair(const NNIOperation &post_nni,
                                       const NNIOperation &pre_nni,
                                       std::optional<size_t> new_tree_id) override;
  // Update the Computation Engine after modifying DAG. Nodes and Edges in reindexer
  // that exceed the prev count gives the location of the new node and edge ids.
  void UpdateEngineAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer) override;
  // Computes scores. Call after Initialize or any Update steps, and before
  // GetTopTreeScores.
  void ComputeScores() override;

  // ** Scoring

  // Get the Top Tree from the DAG with the given edge.
  double GetTopTreeScoreWithEdge(const EdgeId edge_id) const override;
  // Get the Top Tree from the DAG containing the proposed NNI.
  double GetTopTreeScoreWithProposedNNI(
      const NNIOperation &post_nni, const NNIOperation &pre_nni,
      const size_t spare_offset = 0,
      std::optional<BitsetEdgeIdMap> best_edge_map = std::nullopt) override;

  // ** Resize

  // Resize Engine for modified DAG.
  void GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                        std::optional<Reindexer> edge_reindexer) override;
  // Grow engine to handle computing NNIs for all adjacent NNIs.
  // Option to grow engine for computing via reference or via copy. If computing via
  // reference, option whether to use unique temporaries (for testing and computing in
  // parallel).
  void GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                 const bool via_reference = true,
                                 const bool use_unique_temps = true) override;
  // Resize GPEngine to accomodate DAG with given number of nodes and edges.  Option
  // to remap data according to DAG reindexers.  Option to give explicit number of
  // nodes or edges to allocate memory for (this is the only way memory allocation
  // will be decreased).
  void GrowNodeData(const size_t node_count,
                    std::optional<const Reindexer> node_reindexer = std::nullopt,
                    std::optional<const size_t> explicit_alloc = std::nullopt,
                    const bool on_init = false) override;
  void GrowEdgeData(const size_t edge_count,
                    std::optional<const Reindexer> edge_reindexer = std::nullopt,
                    std::optional<const size_t> explicit_alloc = std::nullopt,
                    const bool on_init = false) override;
  // Grow space for storing temporary computation.
  void GrowSpareNodeData(const size_t new_node_spare_count) override;
  void GrowSpareEdgeData(const size_t new_edge_spare_count) override;

  // Copy all edge data from its pre_edge_id to post_edge_id.
  void CopyEdgeData(const EdgeId src_edge_id, const EdgeId dest_edge_id) override;

  // ** Populate PVs

  // Initialize PVs with zeros.
  void ZeroPVs();
  // Populate rootward and leafward PVs.
  void PopulatePVs();
  // Populate P-PVs in a Rootward Pass of DAG.
  void PopulateRootwardPVs();
  // Populate R-PVs in a Leafward Pass of DAG.
  void PopulateLeafwardPVs();

  // ** Branch Length Optimization

  // Initialize branch length handler to prep for optimization.
  // Sets helper functions for optimization methods.
  void InitializeBranchLengthHandler();
  // Perform single pass over DAG of branch length optimization.
  void BranchLengthOptimization(
      std::optional<bool> check_branch_convergence = std::nullopt);
  // Perform branch length optimization on given edge's branch.
  void BranchLengthOptimization(const EdgeId edge_id,
                                const bool check_branch_convergence,
                                const bool update_only = false);

  // Optimization count.
  size_t GetOptimizationCount() { return branch_handler_.GetOptimizationCount(); }
  bool IsFirstOptimization() { return branch_handler_.IsFirstOptimization(); }
  void IncrementOptimizationCount() { branch_handler_.IncrementOptimizationCount(); }
  void ResetOptimizationCount() { branch_handler_.ResetOptimizationCount(); }

  size_t IsOptimizeNewEdges() const { return optimize_new_edges_; }
  void SetOptimizeNewEdges(const bool optimize_new_edges) {
    optimize_new_edges_ = optimize_new_edges;
  }
  size_t GetOptimizationMaxIteration() const { return optimize_max_iter_; }
  void SetOptimizationMaxIteration(const size_t optimize_max_iter) {
    optimize_max_iter_ = optimize_max_iter;
  }

  // ** Access

  PLVEdgeHandler &GetPVs() { return likelihood_pvs_; }
  const PLVEdgeHandler &GetPVs() const { return likelihood_pvs_; }
  EigenMatrixXd &GetMatrix() { return log_likelihoods_; }
  const EigenMatrixXd &GetMatrix() const { return log_likelihoods_; }
  DAGBranchHandler &GetDAGBranchHandler() { return branch_handler_; }
  const DAGBranchHandler &GetDAGBranchHandler() const { return branch_handler_; }

  // ** PV Operations

  struct NNIEdgeIdMap {
    EdgeId central_edge_;
    EdgeId parent_edge_;
    EdgeId sister_edge_;
    EdgeId left_child_edge_;
    EdgeId right_child_edge_;
  };
  struct PrimaryPVIds {
    // For central likelihood.
    PVId parent_rfocal_;
    PVId child_p_;
    // For custom branch lengths.
    PVId child_phatleft_;
    PVId child_phatright_;
    PVId parent_phatsister_;
    PVId parent_rhat_;
    PVId grandparent_rfocal_;
    // For branch length optimization.
    PVId child_rhat_;
    PVId child_rleft_;
    PVId child_rright_;
    PVId parent_phatfocal_;
    PVId parent_rsister_;
    PVId parent_p_;
    PVId grandparent_phatfocal_;
    PVId grandparent_p_;
  };
  struct SecondaryPVIds {
    PVId grandparent_rhat_;
    PVId grandparent_rfocal_;
    PVId grandparent_rsister_;

    PVId parent_rfocal_;
    PVId parent_rsister_;
    PVId child_rleft_;
    PVId child_rright_;

    PVId parent_p_;
    PVId sister_p_;
    PVId leftchild_p_;
    PVId rightchild_p_;

    PVId parent_phatfocal_;
    PVId parent_phatsister_;
    PVId parent_rhat_;
    PVId child_p_;
    PVId child_phatleft_;
    PVId child_phatright_;
    PVId child_rhat_;
  };

  // Get primary PV Ids for corresponding parent/child pair.
  // Gets the P-PV of the child, and the RFocal-PV of the parent.
  // The expected PVs for computing likelihoods.
  std::pair<PVId, PVId> GetPrimaryPVIdsOfEdge(const EdgeId edge_id) const;
  // Get secondary PV Ids for corresponding parent/child pair.
  // Gets the PHatLeft-PV and PHatRight-PV of child, and the the PSister-PV and R-PVs of
  // parent.
  // The expected PVs for computing temp intermediate PVs for proposed NNIs.
  SecondaryPVIds GetSecondaryPVIdsOfEdge(const EdgeId edge_id) const;
  // Remaps secondary PV Ids according to the clade map.
  SecondaryPVIds RemapSecondaryPVIdsForPostNNI(
      const SecondaryPVIds &pre_pvids,
      const NNIOperation::NNICladeArray &clade_map) const;

  // Get temporary PV Ids for intermediate proposed NNI computations.
  PrimaryPVIds GetTempPrimaryPVIdsForProposedNNIs(const size_t spare_offset) const;
  // Get temporary PV Ids for intermediate proposed NNI computations.
  NNIEdgeIdMap GetTempEdgeIdsForProposedNNIs(const size_t spare_offset) const;

  // ** Scoring Helpers

  // Set the P-PVs to match the observed site patterns at the leaves.
  void PopulateLeafPVsWithSitePatterns();
  // Set the R-PVs to the stationary distribution at the root and rootsplits.
  void PopulateRootPVsWithStationaryDistribution();
  // Updates the rootward P-PVs for given node or edge.
  void PopulateRootwardPVForNode(const NodeId node_id);
  void PopulateRootwardPVForEdge(const EdgeId edge_id);
  // Updates the leafward R-PVs for given node or edge.
  void PopulateLeafwardPVForNode(const NodeId node_id);
  void PopulateLeafwardPVForEdge(const EdgeId edge_id);

 protected:
  // ** Scoring Helpers

  // Evolve up the given edge to compute the P-PV of its parent node.
  // Updates the parent's PFocalHat PLV using the child's P PLV.
  void EvolvePPVUpEdge(const EdgeId rootward_edge_id, const EdgeId leafward_edge_id);
  // Evolve down the given edge to compute the R-PV of its child node.
  // Updates the child's RHat PLV using the parent's RFocal PLV.
  void EvolveRPVDownEdge(const EdgeId rootward_edge_id, const EdgeId leafward_edge_id);

  // ** PV Operations
  // Copy data from PV at src_id to dest_id.
  void TakePVValue(const PVId dest_id, const PVId src_id);
  // PV component-wise multiplication of PVs src1 and src2, result stored in dest_id.
  void MultiplyPVs(const PVId dest_id, const PVId src1_id, const PVId src2_id);
  // Compute Likelihood by taking up-to-date parent R-PV and child P-PV.
  void ComputeLikelihood(const EdgeId edge_id, const PVId child_id,
                         const PVId parent_id);
  // Evolve src_id along the branch edge_id and store at dest_id.
  void SetToEvolvedPV(const PVId dest_id, const EdgeId edge_id, const PVId src_id);
  // Evolve src_id along the branch edge_id and multiply with contents of dest_id.
  void MultiplyWithEvolvedPV(const PVId dest_id, const EdgeId edge_id,
                             const PVId src_id);
  // Intermediate computation step for log likelihoods. Stored in temporary variable.
  inline void PreparePerPatternLogLikelihoodsForEdge(const PVId src1_idx,
                                                     const PVId src2_idx) {
    per_pattern_log_likelihoods_ = (GetPVs().GetPV(src1_idx).transpose() *
                                    transition_matrix_ * GetPVs().GetPV(src2_idx))
                                       .diagonal()
                                       .array()
                                       .log();
  }
  // Intermediate computation step for first derivative of log likelihoods. Stored in
  // temporary variable.
  inline void PrepareUnrescaledPerPatternLikelihoodDerivatives(const PVId src1_idx,
                                                               const PVId src2_idx) {
    per_pattern_likelihood_derivatives_ =
        (GetPVs().GetPV(src1_idx).transpose() * derivative_matrix_ *
         GetPVs().GetPV(src2_idx))
            .diagonal()
            .array();
  }
  // Intermediate computation step for second derivative of log likelihoods. Stored in
  // temporary variable.
  inline void PrepareUnrescaledPerPatternLikelihoodSecondDerivatives(
      const PVId src1_idx, const PVId src2_idx) {
    per_pattern_likelihood_second_derivatives_ =
        (GetPVs().GetPV(src1_idx).transpose() * hessian_matrix_ *
         GetPVs().GetPV(src2_idx))
            .diagonal()
            .array();
  }
  // Intermediate computation step for log likelihoods, but without rescaling. Stored in
  // temporary variable.
  inline void PrepareUnrescaledPerPatternLikelihoods(const PVId src1_idx,
                                                     const PVId src2_idx) {
    per_pattern_likelihoods_ = (GetPVs().GetPV(src1_idx).transpose() *
                                transition_matrix_ * GetPVs().GetPV(src2_idx))
                                   .diagonal()
                                   .array();
  }

  // ** Branch Length Optimization Helpers

  // Compute log likelihood and derivative for given edge.
  DoublePair LogLikelihoodAndDerivative(const EdgeId edge_id);
  // Compute log likelihood and first and second derivative for given edge.
  std::tuple<double, double, double> LogLikelihoodAndFirstTwoDerivatives(
      const EdgeId edge_id);
  // Prep temporary transition matrix variable to use given branch length.
  void SetTransitionMatrixToHaveBranchLength(const double branch_length);
  // Prep temporary transition and derivative matrix variables to use given branch
  // length.
  void SetTransitionAndDerivativeMatricesToHaveBranchLength(const double branch_length);
  // Prep temporary transposed transition matrix variable to use given branch length.
  void SetTransitionMatrixToHaveBranchLengthAndTranspose(const double branch_length);

 protected:
  // Partial Vector for storing Likelihood scores.
  PLVEdgeHandler likelihood_pvs_;
  // Tree likelihoods matrix across all sites.
  EigenMatrixXd log_likelihoods_;
  // Branch length parameters for DAG.
  DAGBranchHandler branch_handler_;

  // Determines whether new edges are optimized.
  bool optimize_new_edges_ = true;
  // Number of optimization iterations.
  size_t optimize_max_iter_ = 5;
  // Temporary map of optimized edge lengths.
  std::map<Bitset, double> tmp_optimized_edges;

  // Number of pvs to allocate per node in DAG.
  static constexpr size_t pv_count_per_node_ = PLVTypeEnum::Count;
  // Number of spare nodes needed to be allocated per proposed NNI.
  static constexpr size_t spare_nodes_per_nni_ = 15;
  // Number of spare edges needed to be allocated per proposed NNI.
  static constexpr size_t spare_edges_per_nni_ = 6;

  // ** Substitution Model
  // When we change from JC69Model, check that we are actually doing transpose in
  // leafward calculations.
  JC69Model substitution_model_;
  Eigen::Matrix4d eigenmatrix_ = substitution_model_.GetEigenvectors().reshaped(4, 4);
  Eigen::Matrix4d inverse_eigenmatrix_ =
      substitution_model_.GetInverseEigenvectors().reshaped(4, 4);
  Eigen::Vector4d eigenvalues_ = substitution_model_.GetEigenvalues();
  Eigen::Vector4d stationary_distribution_ = substitution_model_.GetFrequencies();

  // ** Temporaries
  // Stores intermediate computations useful for calculation.
  EigenVectorXd per_pattern_log_likelihoods_;
  EigenVectorXd per_pattern_likelihoods_;
  EigenVectorXd per_pattern_likelihood_derivatives_;
  EigenVectorXd per_pattern_likelihood_derivative_ratios_;
  EigenVectorXd per_pattern_likelihood_second_derivatives_;
  EigenVectorXd per_pattern_likelihood_second_derivative_ratios_;
  Eigen::Vector4d diagonal_vector_;
  Eigen::DiagonalMatrix<double, 4> diagonal_matrix_;
  Eigen::Matrix4d transition_matrix_;
  Eigen::Matrix4d derivative_matrix_;
  Eigen::Matrix4d hessian_matrix_;
};

// TPEngine helper for evaluating Top Trees using Parsimony.
class TPEvalEngineViaParsimony : public TPEvalEngine {
 public:
  TPEvalEngineViaParsimony(TPEngine &tp_engine, const std::string &mmap_path);

  // ** Maintenance

  // Initialize Computation Engine.
  void Initialize() override;
  // Update the Computation Engine after adding Node Pairs to the DAG.
  void UpdateEngineAfterDAGAddNodePair(const NNIOperation &post_nni,
                                       const NNIOperation &pre_nni,
                                       std::optional<size_t> new_tree_id) override;
  // Update the Computation Engine after modifying DAG. Nodes and Edges in reindexer
  // that exceed the prev count gives the location of the new node and edge ids.
  void UpdateEngineAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer) override;
  // Computes scores. Call after Initialize or any Update steps, and before
  // GetTopTreeScores.
  void ComputeScores() override;

  // ** Scoring

  // Get the Top Tree from the DAG with the given edge.
  double GetTopTreeScoreWithEdge(const EdgeId edge_id) const override;
  // Get the Top Tree from the DAG containing the proposed NNI.
  double GetTopTreeScoreWithProposedNNI(
      const NNIOperation &post_nni, const NNIOperation &pre_nni,
      const size_t spare_offset = 0,
      std::optional<BitsetEdgeIdMap> = std::nullopt) override;

  // ** Resize

  // Resize Engine for modified DAG.
  void GrowEngineForDAG(std::optional<Reindexer> node_reindexer,
                        std::optional<Reindexer> edge_reindexer) override;
  // Grow engine to handle computing NNIs for all adjacent NNIs.
  // Option to grow engine for computing via reference or via copy. If computing via
  // reference, option whether to use unique temporaries (for testing and computing in
  // parallel).
  void GrowEngineForAdjacentNNIs(const NNISet &adjacent_nnis,
                                 const bool via_reference = true,
                                 const bool use_unique_temps = true) override;
  // Resize GPEngine to accomodate DAG with given number of nodes and edges.  Option
  // to remap data according to DAG reindexers.  Option to give explicit number of
  // nodes or edges to allocate memory for (this is the only way memory allocation
  // will be decreased).
  void GrowNodeData(const size_t node_count,
                    std::optional<const Reindexer> node_reindexer = std::nullopt,
                    std::optional<const size_t> explicit_alloc = std::nullopt,
                    const bool on_init = false) override;
  void GrowEdgeData(const size_t edge_count,
                    std::optional<const Reindexer> edge_reindexer = std::nullopt,
                    std::optional<const size_t> explicit_alloc = std::nullopt,
                    const bool on_init = false) override;
  // Grow space for storing temporary computation.
  void GrowSpareNodeData(const size_t new_node_spare_count) override;
  void GrowSpareEdgeData(const size_t new_edge_spare_count) override;

  // Copy all edge data from its pre_edge_id to post_edge_id.
  void CopyEdgeData(const EdgeId src_edge_id, const EdgeId dest_edge_id) override;

  // ** Populate PVs

  // Initialize PVs with zero.
  void ZeroPVs();
  // Populate rootward and leafward PVs.
  void PopulatePVs();
  // Populate P-PVs in a Rootward Pass of DAG.
  void PopulateRootwardPVs();
  // Populate R-PVs in a Leafward Pass of DAG.
  void PopulateLeafwardPVs();

  // ** Access

  PSVEdgeHandler &GetPVs() { return parsimony_pvs_; }
  const PSVEdgeHandler &GetPVs() const { return parsimony_pvs_; }

 protected:
  // ** Scoring Helpers

  // Compute the rootward P-PVs for given node or edge.
  void PopulateRootwardParsimonyPVForNode(const NodeId node_id);
  void PopulateRootwardParsimonyPVForEdge(const EdgeId edge_id);
  // Compute the leafward R-PVs for given node or edge.
  void PopulateLeafwardParsimonyPVForNode(const NodeId node_id);
  void PopulateLeafwardParsimonyPVForEdge(const EdgeId edge_id);
  // Set the P-PVs to match the observed site patterns at the leaves.
  void PopulateLeafParsimonyPVsWithSitePatterns();
  // Calculate the PV for a given parent-child pair.
  EigenVectorXd ParentPartial(EigenVectorXd child_partials);
  // Sum P-PVs for right and left children of node 'node_id'
  // In this case, we get the full P-PVs of the given node after all P-PVs
  // have been concatenated into one SankoffPartialVector.
  EigenVectorXd TotalPPartial(const EdgeId edge_id, const size_t site_idx);
  EigenVectorXd TotalPPartial(const PVId edge_pleft_pvid, const PVId edge_pright_pvid,
                              const size_t site_idx);
  // Populate rootward P-PVs for given edge.
  // Updates parent's Pleft PV with sum of left child P PVs and parent's PRight PV with
  // sum of right child P PVs.
  void PopulateRootwardParsimonyPVForEdge(const EdgeId parent_id,
                                          const EdgeId left_child_id,
                                          const EdgeId right_child_id);
  // Populate leafward Q-PVs for given edge.
  // Updates parent's Q PV by combining with sister's sum of P PVs.
  void PopulateLeafwardParsimonyPVForEdge(const EdgeId parent_id,
                                          const EdgeId left_child_id,
                                          const EdgeId right_child_id);
  // Calculates parsimony score on given edge.
  // Takes the minimum element after taking sum of edge's P and Q PVs.
  double ParsimonyScore(const EdgeId edge_id);
  double ParsimonyScore(const PVId edge_q_pvid, const PVId edge_pleft_pvid,
                        const PVId edge_pright_pvid);

 protected:
  // Partial Vector for computing Parsimony scores.
  PSVEdgeHandler parsimony_pvs_;
  static constexpr size_t pv_count_per_node_ = PSVTypeEnum::Count;
  // Stores intermediate computations.
  // Internal "temporaries" useful for calculation.
  SankoffMatrix parsimony_cost_matrix_;

  static constexpr size_t state_count_ = 4;
  static constexpr double big_double_ = static_cast<double>(INT_MAX);

  // Number of spare nodes needed to be allocated per proposed NNI.
  static constexpr size_t spare_nodes_per_nni_ = 2;
  // Number of spare edges needed to be allocated per proposed NNI.
  static constexpr size_t spare_edges_per_nni_ = 5;
};
