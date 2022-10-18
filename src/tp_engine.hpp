// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// TPEngine runs the Top-Pruning method for systematic exploration. This method
// procedurally explores NNIs adjacent to the DAG in tree space.  NNIs are then
// evaluated by finding the best scoring tree contained in the DAG which has a
// given edge.  Scoring of trees can be done through likelihoods or parsimony.

#pragma once

#include "sugar.hpp"
#include "gp_dag.hpp"
#include "graft_dag.hpp"
#include "pv_handler.hpp"
#include "choice_map.hpp"
#include "nni_operation.hpp"
#include "sankoff_handler.hpp"
#include "dag_data.hpp"
#include "optimization.hpp"

class TPEngine {
 public:
  TPEngine(GPDAG &dag, SitePattern &site_pattern,
           const std::string &mmap_likelihood_path, bool using_likelihood,
           const std::string &mmap_parsimony_path, bool using_parsimony,
           bool use_gradients = true);

  // ** Access

  EigenVectorXd &GetTopTreeLikelihoodsPerEdge() {
    return top_tree_log_likelihoods_per_edge_;
  }
  EigenVectorXd &GetTopTreeParsimonyPerEdge() { return top_tree_parsimony_per_edge_; }

  // ** General Scoring

  // Intialize choice map naively by setting first encountered edge for each
  void InitializeChoiceMap();
  // Get top-scoring tree in DAG containing given edge.
  Node::NodePtr GetTopTreeTopologyWithEdge(const EdgeId edge_id) const;
  // Get score of top-scoring tree in DAG containing given edge.
  double GetTopTreeScore(const EdgeId edge_id) const;

  // ** Scoring by Likelihood

  // Initialize ChoiceMap and populate PVs for entire DAG by Likelihood.
  void InitializeLikelihood();
  // Fetch likelihood of top tree with given edge.  Assumed likelihoods have already
  // been computed.
  double GetTopTreeLikelihoodWithEdge(const EdgeId edge_id);
  // Compute top tree likelihoods for all edges in DAG. Result stored in
  // log_likelihoods_per_edge_ vector..
  void ComputeLikelihoods();
  // After adding an NNI to the DAG, update the out-of-date likelihoods.
  void UpdateLikelihoodsAfterDAGAddNodePair(
      const NNIOperation &post_nni, const NNIOperation &pre_nni,
      std::optional<size_t> new_tree_id = std::nullopt);
  // Find the likelihood of a proposed NNI (an NNI not currently stored in the DAG but
  // adjacent to nodes in the current DAG).  Computations are done in place, leveraging
  // pre-existing PVs.
  double GetTopTreeLikelihoodWithProposedNNI(const NNIOperation &post_nni,
                                             const NNIOperation &pre_nni,
                                             const size_t spare_offset = 0);

  // ** Scoring by Parsimony

  // Initialize ChoiceMap and populate PVs for entire DAG by Parsimony.
  void InitializeParsimony();
  // Fetch parsimony of top tree with given edge.  Assumed parsimony have already
  // been computed.
  double GetTopTreeParsimonyWithEdge(const EdgeId edge_id);
  // Compute top tree parsimony for all edges in DAG. Result stored in
  // log_parsimony_per_edge_ vector.
  void ComputeParsimonies();
  // After adding an NNI to the DAG, update the out-of-date likelihoods.
  void UpdateParsimoniesAfterDAGAddNodePair(
      const NNIOperation &post_nni, const NNIOperation &pre_nni,
      std::optional<size_t> new_tree_id = std::nullopt);

  // ** Parameter Data

  // Resize GPEngine to accomodate DAG with given number of nodes and edges.  Option
  // to remap data according to DAG reindexers.  Option to give explicit number of
  // nodes or edges to allocate memory for (this is the only way memory allocation
  // will be decreased).
  void GrowNodeData(const size_t node_count,
                    std::optional<const Reindexer> node_reindexer = std::nullopt,
                    std::optional<const size_t> explicit_allocation = std::nullopt,
                    const bool on_initialization = false);
  void GrowEdgeData(const size_t edge_count,
                    std::optional<const Reindexer> edge_reindexer = std::nullopt,
                    std::optional<const size_t> explicit_allocation = std::nullopt,
                    const bool on_intialization = false);
  // Remap node and edge-based data according to reordering of DAG nodes and edges.
  void ReindexNodeData(const Reindexer &node_reindexer, const size_t old_node_count);
  void ReindexEdgeData(const Reindexer &edge_reindexer, const size_t old_edge_count);
  // Grow space for storing temporary computation.
  void GrowSpareNodeData(const size_t new_node_spare_count);
  void GrowSpareEdgeData(const size_t new_edge_spare_count);

  // ** Tree Collection

  // Set choice map by taking the first occurrence of each PCSP edge from collection.
  void SetChoiceMapByTakingFirst(const RootedTreeCollection &tree_collection,
                                 const BitsetSizeMap &edge_indexer);
  // Set branch lengths by taking the first occurrance of each PCSP edge from
  // collection.
  void SetBranchLengthByTakingFirst(const RootedTreeCollection &tree_collection,
                                    const BitsetSizeMap &edge_indexer);

  // ** Counts

  // Node Counts
  size_t GetNodeCount() const { return node_count_; };
  size_t GetSpareNodeCount() const { return node_spare_count_; }
  size_t GetAllocatedNodeCount() const { return node_alloc_; }
  size_t GetPaddedNodeCount() const { return GetNodeCount() + GetSpareNodeCount(); };
  void SetNodeCount(const size_t node_count) {
    node_count_ = node_count;
    if (using_likelihoods_) {
      likelihood_pvs_.SetCount(node_count);
    }
    if (using_parsimony_) {
      parsimony_pvs_.SetCount(node_count);
    }
  }
  void SetSpareNodeCount(const size_t node_spare_count) {
    node_spare_count_ = node_spare_count;
    if (using_likelihoods_) {
      likelihood_pvs_.SetSpareCount(node_spare_count);
    }
    if (using_parsimony_) {
      parsimony_pvs_.SetSpareCount(node_spare_count);
    }
  }
  void SetAllocatedNodeCount(const size_t node_alloc) {
    node_alloc_ = node_alloc;
    if (using_likelihoods_) {
      likelihood_pvs_.SetAllocatedCount(node_alloc);
    }
    if (using_parsimony_) {
      parsimony_pvs_.SetAllocatedCount(node_alloc);
    }
  }

  // Edge Counts
  size_t GetEdgeCount() const { return edge_count_; };
  size_t GetSpareEdgeCount() const { return edge_spare_count_; };
  size_t GetAllocatedEdgeCount() const { return edge_alloc_; };
  size_t GetPaddedEdgeCount() const { return GetEdgeCount() + GetSpareEdgeCount(); };
  size_t GetSpareEdgeIndex(const size_t edge_offset) const {
    const size_t edge_scratch_size = GetPaddedEdgeCount() - GetEdgeCount();
    Assert(edge_offset < edge_scratch_size,
           "Requested edge_offset outside of allocated scratch space.");
    return edge_offset + GetEdgeCount();
  }
  void SetEdgeCount(const size_t edge_count) { edge_count_ = edge_count; }
  void SetSpareEdgeCount(const size_t edge_spare_count) {
    edge_spare_count_ = edge_spare_count;
  }
  void SetAllocatedEdgeCount(const size_t edge_alloc) { edge_alloc_ = edge_alloc; }

  size_t GetInputTreeCount() const { return input_tree_count_; }

  // ** Access

  EigenConstMatrixXdRef GetLikelihoodMatrix() {
    return log_likelihoods_.block(0, 0, GetNodeCount(), log_likelihoods_.cols());
  }
  PLVEdgeHandler &GetLikelihoodPVs() { return likelihood_pvs_; }
  PSVEdgeHandler &GetParsimonyPVs() { return parsimony_pvs_; }
  EigenVectorXd &GetBranchLengths() { return branch_lengths_; }
  std::vector<size_t> &GetTreeSource() { return tree_source_; }
  EigenVectorXd &GetTopTreeLikelihoods() { return top_tree_log_likelihoods_per_edge_; }
  EigenVectorXd &GetTopTreeParsimonies() { return top_tree_parsimony_per_edge_; }

  void SetBranchLengths(EigenVectorXd branch_lengths) {
    Assert(size_t(branch_lengths.size()) == dag_.EdgeCountWithLeafSubsplits(),
           "Size mismatch in GPEngine::SetBranchLengths.");
    branch_lengths_.segment(0, dag_.EdgeCountWithLeafSubsplits()) = branch_lengths;
  }

  // ** I/O

  std::string LikelihoodPVToString(const PVId pv_id) const;
  std::string LogLikelihoodMatrixToString() const;
  std::string ParsimonyPVToString(const PVId pv_id) const;

 protected:
  // ** Scoring

  // Update edge and node data by copying over from pre-NNI to post-NNI.
  void CopyOverEdgeDataFromPreNNIToPostNNI(
      const NNIOperation &post_nni, const NNIOperation &pre_nni,
      std::optional<size_t> new_tree_id = std::nullopt);
  // Find the edge from the best scoring tree that is adjacent to given node in
  // the given direction.
  // Accomplished by iterating over all adjacent edges using tree_source_ edge
  // map, which gives the best tree id using a given edge. The best tree is
  // expected to be the earliest found in the tree collection, aka smallest
  // tree id. The adjacent edge that comes from the best tree is chosen.
  EdgeId FindBestEdgeAdjacentToNode(const NodeId node_id,
                                    const Direction direction) const;
  EdgeId FindBestEdgeAdjacentToNode(const NodeId node_id, const Direction direction,
                                    const SubsplitClade clade) const;

  // ** Scoring by Likelihoods

  // Compute the rootward P-PVs for given node.
  void PopulateRootwardLikelihoodPVForNode(const NodeId node_id);
  // Compute the leafward R-PVs for given node.
  void PopulateLeafwardLikelihoodPVForNode(const NodeId node_id);
  // Set the P-PVs to match the observed site patterns at the leaves.
  void PopulateLeafLikelihoodPVsWithSitePatterns();
  // Set the R-PVs to the stationary distribution at the root and rootsplits.
  void PopulateRootLikelihoodPVsWithStationaryDistribution();
  // Evolve up the given edge to compute the P-PV of its parent node.
  void EvolveLikelihoodPPVUpEdge(const EdgeId parent_edge_id,
                                 const EdgeId child_edge_id);
  // Evolve down the given edge to compute the R-PV of its child node.
  void EvolveLikelihoodRPVDownEdge(const EdgeId parent_edge_id,
                                   const EdgeId child_edge_id);

  // ** Scoring by Parsimony

  // Compute the rootward P-PVs for given node.
  void PopulateRootwardParsimonyPVForNode(const NodeId node_id);
  // Compute the leafward R-PVs for given node.
  void PopulateLeafwardParsimonyPVForNode(const NodeId node_id);
  // Set the P-PVs to match the observed site patterns at the leaves.
  void PopulateLeafParsimonyPVsWithSitePatterns();
  // Calculate the PV for a given parent-child pair.
  EigenVectorXd ParentPartial(EigenVectorXd child_partials);
  // Sum P-PVs for right and left children of node 'node_id'
  // In this case, we get the full P-PVs of the given node after all P-PVs
  // have been concatenated into one SankoffPartialVector.
  EigenVectorXd TotalPPartial(EdgeId edge_id, size_t site_idx);
  // Populate rootward R-PVs for given edge.
  void PopulateRootwardParsimonyPVForEdge(const EdgeId parent_id,
                                          const EdgeId left_child_id,
                                          const EdgeId right_child_id);
  // Populate leafward P-PVs for given edge.
  void PopulateLeafwardParsimonyPVForEdge(const EdgeId parent_id,
                                          const EdgeId left_child_id,
                                          const EdgeId right_child_id);
  // Calculates parsimony score on given edge across all sites.
  double ParsimonyScore(const EdgeId edge_id);

  // ** PV Operations for Likelihoods

  // Assign PV at src_id to dest_id.
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
  // Prepare to evolve along given branch length. Stored in temporary variables
  // diagonal_matrix_ and transition_matrix_.
  void SetTransitionMatrixToHaveBranchLength(const double branch_length);
  // Intermediate likelihood computation step. Stored in temporary variable
  // per_pattern_log_likelihoods_.
  void PreparePerPatternLogLikelihoodsForEdge(const PVId src1_id, const PVId src2_id);

  // ** DAG
  // Un-owned reference DAG.
  GPDAG &dag_;
  // For adding temporary NNIs to DAG.
  std::unique_ptr<GraftDAG> graft_dag_;

  // ** Parameters
  // Total number of nodes in DAG. Determines sizes of data vectors indexed on nodes.
  size_t node_count_ = 0;
  size_t node_alloc_ = 0;
  size_t node_spare_count_ = 2;
  // Total number of edges in DAG. Determines sizes of data vectors indexed on edges.
  size_t edge_count_ = 0;
  size_t edge_alloc_ = 0;
  size_t edge_spare_count_ = 3;
  // Growth factor when reallocating data.
  constexpr static double resizing_factor_ = 2.0;

  // Tree likelihoods matrix across all sites.
  EigenMatrixXd log_likelihoods_;
  // Top tree log likelihood per edge.
  EigenVectorXd top_tree_log_likelihoods_per_edge_;

  // Tree parsimonies across all sites.
  EigenMatrixXd parsimonies_;
  // Top tree parsimony per edge.
  EigenVectorXd top_tree_parsimony_per_edge_;

  // ** Branch length parameters for DAG.
  EigenVectorXd branch_lengths_;
  EigenVectorXd branch_length_differences_;
  // Initial branch length during first branch length opimization.
  static constexpr double default_branch_length_ = 0.1;
  // Absolute lower bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double min_log_branch_length_ = -13.9;
  // Absolute upper bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double max_log_branch_length_ = 1.1;
  // Method used for performing optimization.
  Optimization::OptimizationMethod optimization_method_;
  // Precision used for checking convergence of branch length optimization.
  // In the non-Brent optimization methods, significant digits will be used to
  // determine convergence through relative tolerance, i.e. measuring difference
  // from previous branch length values until the absolute difference is below
  // 10^(-significant_digits_for_optimization_).
  // Brent optimization does not define convergence through relative tolerance,
  // rather convergence based on tightness of the brackets that it adapts during
  // optimization. This variable thus represents the "number of bits precision to which
  // the minima should be found". When testing on sample datasets, we found that setting
  // the value to 10 was a good compromise between speed and precision for Brent.
  // See more on Brent optimization here:
  // https://www.boost.org/doc/libs/1_79_0/libs/math/doc/html/math_toolkit/brent_minima.html
  int significant_digits_for_optimization_ = 10;
  double relative_tolerance_for_optimization_ = 1e-4;
  double denominator_tolerance_for_newton_ = 1e-10;
  double step_size_for_optimization_ = 5e-4;
  double step_size_for_log_space_optimization_ = 1.0005;
  // Number of iterations allowed for branch length optimization.
  size_t max_iter_for_optimization_ = 1000;
  double branch_length_difference_threshold_ = 1e-15;

  // Observed leave states.
  SitePattern site_pattern_;

  // Tree id where branch_length and choice_map is sourced.
  // TreeCollection is expected to be ordered from highest to lowest scoring, so lower
  // tree id means better scoring tree.
  std::vector<size_t> tree_source_;
  // Total number of trees used to construct the DAG.
  size_t input_tree_count_ = 0;

  // ** Scoring
  // Partial Vector for storing Likelihood scores.
  PLVEdgeHandler likelihood_pvs_;
  bool using_likelihoods_;
  // Partial Vector for storing Parsimony scores.
  PSVEdgeHandler parsimony_pvs_;
  SankoffMatrix parsimony_cost_matrix_;
  bool using_parsimony_;
  // ChoiceMap for find top-scoring tree containing any given branch.
  ChoiceMap choice_map_;

  // ** Temporaries
  // Stores intermediate values for computation.
  EigenVectorXd per_pattern_log_likelihoods_;

  // ** Substitution Model
  // When we change from JC69Model, check that we are actually doing transpose in
  // leafward calculations.
  JC69Model substitution_model_;
  Eigen::Matrix4d eigenmatrix_ = substitution_model_.GetEigenvectors().reshaped(4, 4);
  Eigen::Matrix4d inverse_eigenmatrix_ =
      substitution_model_.GetInverseEigenvectors().reshaped(4, 4);
  Eigen::Vector4d eigenvalues_ = substitution_model_.GetEigenvalues();
  Eigen::Vector4d diagonal_vector_;
  Eigen::DiagonalMatrix<double, 4> diagonal_matrix_;
  Eigen::Matrix4d transition_matrix_;
  Eigen::Vector4d stationary_distribution_ = substitution_model_.GetFrequencies();
  EigenVectorXd site_pattern_weights_;

  static constexpr size_t state_count_ = 4;
  static constexpr double big_double_ = static_cast<double>(INT_MAX);
};
