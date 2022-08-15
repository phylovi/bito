// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// TPEngine runs the Top-Pruning method for systematic exploration. This method
// procedurally explores NNIs adjacent to the DAG in tree space.  NNIs are then
// evaluated by finding the best scoring tree contained in the DAG which has a
// given edge.  Scoring of trees can be done through likelihoods or parsimony.

#include "sugar.hpp"
#include "gp_dag.hpp"
#include "graft_dag.hpp"
#include "pv_handler.hpp"
#include "choice_map.hpp"
#include "nni_operation.hpp"

using PVId = size_t;

class TPEngine {
 public:
  TPEngine(GPDAG &dag, SitePattern &site_pattern, const std::string &mmap_file_path,
           bool using_likelihood, bool using_parsimony);

  // ** General Scoring

  // Intialize choice map naively by setting first encountered edge for each
  void InitializeChoiceMap();
  // Get top-scoring tree in DAG containing given edge.
  Node::NodePtr GetTopTreeTopologyWithEdge(const EdgeId edge_id) const;
  // Get score of top-scoring tree in DAG containing given edge.
  double GetTopTreeScore(const EdgeId edge_id) const;

  // ** Scoring by Likelihood

  // Initialize ChoiceMap for entire DAG by Likelihood.
  void InitializeLikelihood();
  // Update ChoiceMap with NNI.
  void UpdateDAGAfterAddNodePairByLikelihood(const NNIOperation &nni_op);
  // Fetch likelihood of top tree with given edge.  Assumed likelihoods have already
  // been computed.
  double GetTopTreeLikelihoodWithEdge(const EdgeId edge_id);
  // Compute top tree likelihoods for all edges in DAG. Result stored in
  // log_likelihoods_ matrix.
  void ComputeLikelihoods();

  // ** Scoring by Parsimony

  // Initialize ChoiceMap for entire DAG by Parsimony.
  void InitializeParsimony();
  // Update ChoiceMap with NNI.
  void UpdateDAGAfterAddNodePairByParsimony(const NNIOperation &nni_op);

  // ** Parameters

  // Resize GPEngine to accomodate DAG with given number of nodes and edges.  Option to
  // remap data according to DAG reindexers.  Option to give explicit number of nodes or
  // edges to allocate memory for (this is the only way memory allocation will be
  // decreased).
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

  // Apply function to edges descending from each node on each rooted tree for all trees
  // in collection.
  using FunctionOnTreeNodeByEdge = std::function<void(
      const EdgeId, const Bitset &, const RootedTree &, const size_t, const Node *)>;
  void FunctionOverRootedTreeCollection(
      FunctionOnTreeNodeByEdge function_on_tree_node_by_edge,
      const RootedTreeCollection &tree_collection, const BitsetSizeMap &indexer);

  // Set choice map by taking the first occurrance of each PCSP edge from collection.
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
      likelihood_pvs_.SetNodeCount(node_count);
    }
    if (using_parsimony_) {
      parsimony_pvs_.SetNodeCount(node_count);
    }
  }
  void SetSpareNodeCount(const size_t node_spare_count) {
    node_spare_count_ = node_spare_count;
    if (using_likelihoods_) {
      likelihood_pvs_.SetSpareNodeCount(node_spare_count);
    }
    if (using_parsimony_) {
      parsimony_pvs_.SetSpareNodeCount(node_spare_count);
    }
  }
  void SetAllocatedNodeCount(const size_t node_alloc) {
    node_alloc_ = node_alloc;
    if (using_likelihoods_) {
      likelihood_pvs_.SetAllocatedNodeCount(node_alloc);
    }
    if (using_parsimony_) {
      parsimony_pvs_.SetAllocatedNodeCount(node_alloc);
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
  PLVHandler &GetLikelihoodPVs() { return likelihood_pvs_; }
  PSVHandler &GetParsimonyPVs() { return parsimony_pvs_; }
  EigenVectorXd &GetBranchLengths() { return branch_lengths_; }
  std::vector<size_t> &GetTreeSource() { return tree_source_; }

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
  // ** Scoring by Likelihoods

  void PopulateRootwardLikelihoodPVForNode(const NodeId node_id);
  void PopulateLeafwardLikelihoodPVForNode(const NodeId node_id);
  void PopulateLeafLikelihoodPVsWithSitePatterns();
  void PopulateRootLikelihoodPVsWithStationaryDistribution();
  EdgeId BestEdgeAdjacentToNode(const NodeId node_id, const Direction direction) const;
  void EvolveLikelihoodPPVUpEdge(const EdgeId edge_id);
  void EvolveLikelihoodRPVDownEdge(const EdgeId edge_id);

  // ** Scoring by Parsimony

  void PopulateLeafParsimonyPVsWithSitePatterns();

  // ** Partial Vector Operations

  // Assign PV at src_id to dest_id.
  void TakePVValue(const PVId dest_id, const PVId src_id);
  // PV component-wise multiplication of PVs src1 and src2, result stored in dest_id.
  void MultiplyPVs(const PVId dest_id, const PVId src1_id, const PVId src2_id);
  // Compute Likelihood by taking up-to-date parent R-PLV and child P-PLV.
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
  // Total number of trees used to construct the DAG.
  size_t input_tree_count_ = 0;
  // Growth factor when reallocating data.
  constexpr static double resizing_factor_ = 2.0;

  // Tree likelihoods matrix across all sites.
  EigenMatrixXd log_likelihoods_;
  // Top tree log likelihood per edge.
  EigenVectorXd top_tree_log_likelihoods_per_edge_;

  // Branch length parameters for DAG.
  EigenVectorXd branch_lengths_;
  // Initial branch length during first branch length opimization.
  static constexpr double default_branch_length_ = 0.1;
  // Absolute lower bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double min_log_branch_length_ = -13.9;
  // Absolute upper bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double max_log_branch_length_ = 1.1;

  // Observed leave states.
  SitePattern site_pattern_;

  // Tree id where branch_length and choice_map is sourced.
  // TreeCollection is expected to be ordered from highest to lowest scoring, so lower
  // ID means higher priority tree.
  std::vector<size_t> tree_source_;

  // ** Scoring
  // Partial Vector for storing Likelihood scores.
  PLVHandler likelihood_pvs_;
  bool using_likelihoods_;
  // Partial Vector for storing Parsimony scores.
  PSVHandler parsimony_pvs_;
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
};
