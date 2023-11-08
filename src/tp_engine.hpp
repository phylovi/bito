// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// TPEngine uses the Top-Pruning method for evaluating the edges of the DAG. Each edge
// is evaluated according to its representative "Top Tree", that is, the best scoring
// tree of all possible trees contained within the DAG that include that edge.  The
// topology of each Top Tree is stored per edge in the choice map. Tree/edge scoring is
// facilitated by the helper TPEvalEngine class.

#pragma once

#include "sugar.hpp"
#include "gp_dag.hpp"
#include "graft_dag.hpp"
#include "pv_handler.hpp"
#include "tp_choice_map.hpp"
#include "nni_operation.hpp"
#include "dag_branch_handler.hpp"
#include "optimization.hpp"
#include "substitution_model.hpp"
#include "tp_evaluation_engine.hpp"

enum class TPEvalEngineType { LikelihoodEvalEngine, ParsimonyEvalEngine };
static const inline size_t TPEvalEngineTypeCount = 2;
class TPEvalEngineTypeEnum
    : public EnumWrapper<TPEvalEngineType, size_t, TPEvalEngineTypeCount,
                         TPEvalEngineType::LikelihoodEvalEngine,
                         TPEvalEngineType::ParsimonyEvalEngine> {
 public:
  static inline const std::string Prefix = "TPEvalEngineType";
  static inline const Array<std::string> Labels = {
      {"LikelihoodEvalEngine", "ParsimonyEvalEngine"}};

  static std::string ToString(const TPEvalEngineType e) {
    std::stringstream ss;
    ss << Prefix << "::" << Labels[e];
    return ss.str();
  }
  friend std::ostream &operator<<(std::ostream &os, const TPEvalEngineType e) {
    os << ToString(e);
    return os;
  }
};

using BitsetEdgeIdMap = std::unordered_map<Bitset, EdgeId>;
using BitsetBitsetMap = std::unordered_map<Bitset, Bitset>;
using NNINNIMap = std::unordered_map<NNIOperation, NNIOperation>;
using NNIAdjNodeIdMap = NNIAdjacentMap<std::pair<NodeId, NodeId>>;
using NNIAdjBitsetEdgeIdMap = NNIAdjacentMap<std::pair<Bitset, EdgeId>>;

class TPEngine {
 public:
  TPEngine(GPDAG &dag, SitePattern &site_pattern);
  TPEngine(GPDAG &dag, SitePattern &site_pattern,
           std::optional<std::string> mmap_likelihood_path,
           std::optional<std::string> mmap_parsimony_path,
           std::optional<const RootedTreeCollection> tree_collection = std::nullopt,
           std::optional<const BitsetSizeMap> edge_indexer = std::nullopt);

  // ** Comparators

  // Compares DAG and EdgeChoices. Note: only tests for equality, not a well-ordered
  // comparator.
  static int Compare(const TPEngine &lhs, const TPEngine &rhs,
                     const bool is_quiet = true);
  friend bool operator==(const TPEngine &lhs, const TPEngine &rhs);

  // ** Access

  GPDAG &GetDAG() { return *dag_; }
  const GPDAG &GetDAG() const { return *dag_; }
  GraftDAG &GetGraftDAG() { return *graft_dag_; }
  const GraftDAG &GetGraftDAG() const { return *graft_dag_; }
  SitePattern &GetSitePattern() { return site_pattern_; }
  const SitePattern &GetSitePattern() const { return site_pattern_; }
  EigenVectorXd &GetSitePatternWeights() { return site_pattern_weights_; }
  const EigenVectorXd &GetSitePatternWeights() const { return site_pattern_weights_; }

  TPChoiceMap &GetChoiceMap() { return choice_map_; }
  const TPChoiceMap &GetChoiceMap() const { return choice_map_; }
  TPChoiceMap::EdgeChoice &GetChoiceMap(const EdgeId edge_id) {
    return GetChoiceMap().GetEdgeChoice(edge_id);
  }
  const TPChoiceMap::EdgeChoice &GetChoiceMap(const EdgeId edge_id) const {
    return GetChoiceMap().GetEdgeChoice(edge_id);
  }

  std::vector<TreeId> &GetTreeSource() { return tree_source_; }
  const std::vector<TreeId> &GetTreeSource() const { return tree_source_; }
  void SetTreeSource(const std::vector<TreeId> tree_source) {
    tree_source_ = tree_source;
  }
  TreeId &GetTreeSource(const EdgeId edge_id) {
    return GetTreeSource()[edge_id.value_];
  }
  const TreeId &GetTreeSource(const EdgeId edge_id) const {
    return GetTreeSource()[edge_id.value_];
  }

  // ** Counts

  // Node Counts
  size_t GetNodeCount() const { return node_count_; };
  size_t GetSpareNodeCount() const { return node_spare_count_; }
  size_t GetAllocatedNodeCount() const { return node_alloc_; }
  size_t GetPaddedNodeCount() const { return GetNodeCount() + GetSpareNodeCount(); };
  void SetNodeCount(const size_t node_count) { node_count_ = node_count; }
  void SetSpareNodeCount(const size_t node_spare_count) {
    node_spare_count_ = node_spare_count;
  }
  void SetAllocatedNodeCount(const size_t node_alloc) { node_alloc_ = node_alloc; }

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

  size_t GetSpareNodesPerNNI() const { return spare_nodes_per_nni_; }
  size_t GetSpareEdgesPerNNI() const { return spare_edges_per_nni_; }

  double GetResizingFactor() const { return resizing_factor_; }
  size_t GetInputTreeCount() const { return input_tree_count_; }
  TreeId GetMaxTreeId() const { return TreeId(tree_counter_); }
  TreeId GetNextTreeId() const { return TreeId(GetInputTreeCount()); }

  // ** Maintenance

  // Initialize engine parts: choice map, eval engine, etc.
  void Initialize();
  // Update engine after updating the DAG.
  void UpdateAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer,
      bool is_quiet = true);
  // Resize GPEngine to accomodate DAG with given number of nodes and edges.  Option
  // to remap data according to DAG reindexers.  Option to give explicit number of
  // nodes or edges to allocate memory for (this is the only way memory allocation
  // will be decreased).
  void GrowNodeData(const size_t node_count,
                    std::optional<const Reindexer> node_reindexer = std::nullopt,
                    std::optional<const size_t> explicit_alloc = std::nullopt,
                    const bool on_init = false);
  void GrowEdgeData(const size_t edge_count,
                    std::optional<const Reindexer> edge_reindexer = std::nullopt,
                    std::optional<const size_t> explicit_alloc = std::nullopt,
                    const bool on_intialization = false);
  // Remap node and edge-based data according to reordering of DAG nodes and edges.
  void ReindexNodeData(const Reindexer &node_reindexer, const size_t old_node_count);
  void ReindexEdgeData(const Reindexer &edge_reindexer, const size_t old_edge_count);
  // Grow space for storing temporary computation.
  void GrowSpareNodeData(const size_t new_node_spare_count);
  void GrowSpareEdgeData(const size_t new_edge_spare_count);

  // Update edge and node data by copying over from pre-NNI to post-NNI.
  using CopyEdgeDataFunc = std::function<void(const EdgeId, const EdgeId)>;
  void CopyOverEdgeDataFromPreNNIToPostNNI(
      const NNIOperation &post_nni, const NNIOperation &pre_nni,
      CopyEdgeDataFunc copy_data_func,
      std::optional<size_t> new_tree_id = std::nullopt);

  // ** Choice Map

  // Intialize choice map naively by setting first encountered edge for each
  void InitializeChoiceMap();
  // Update choice map after modifying DAG.
  void UpdateChoiceMapAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer);

  // Set tree source for each edge, by taking the first occurrence of each PCSP edge
  // from input trees.
  void SetTreeSourceByTakingFirst(const RootedTreeCollection &tree_collection,
                                  const BitsetSizeMap &edge_indexer);
  // Set each edge's choice map, by either:
  // True: the PCSP heuristic, False: the Subsplit heuristic.
  void SetChoiceMapByTakingFirst(const RootedTreeCollection &tree_collection,
                                 const BitsetSizeMap &edge_indexer,
                                 const bool use_subsplit_heuristic = true);

  // Update an individual edge's choice map using the tree source. Naively takes first
  // adjacent edge.
  void UpdateEdgeChoiceByTakingFirstTree(const EdgeId edge_id);
  // Update an individual edge's choice map using the tree source. Trees added to the
  // DAG first recieve highest priority.
  void UpdateEdgeChoiceByTakingHighestPriorityTree(const EdgeId edge_id);
  // Update an individual edge's choice map using the tree source. Examines all
  // available edge combinations an takes adjacent edge that results in max score.
  void UpdateEdgeChoiceByTakingHighestScoringTree(const EdgeId edge_id);

  // ** Proposed NNIs

  // Get highest priority pre-NNI in DAG for given post-NNI.
  NNIOperation FindHighestPriorityNeighborNNIInDAG(const NNIOperation &nni) const;
  std::tuple<EdgeId, EdgeId, EdgeId> FindHighestPriorityAdjacentNodeId(
      const NodeId node_id) const;
  // Builds a map of adjacent edges from pre-NNI to post-NNI.
  std::unordered_map<EdgeId, EdgeId> BuildAdjacentEdgeMapFromPostNNIToPreNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni) const;

  // Map edge ids in pre-NNI edge choice according to the swap in post-NNI clade map.
  TPChoiceMap::EdgeChoice RemapEdgeChoiceFromPreNNIToPostNNI(
      const TPChoiceMap::EdgeChoice &pre_choice,
      const NNIOperation::NNICladeArray &clade_map) const;
  // Create remapped edge choices from pre-NNI to post-NNI.
  TPChoiceMap::EdgeChoice GetRemappedEdgeChoiceFromPreNNIToPostNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni) const;

  // Get the average of the edges below parent and above child.
  double GetAvgLengthOfAdjEdges(
      const NodeId parent_node_id, const NodeId child_node_id,
      const std::optional<size_t> prev_node_count = std::nullopt,
      const std::optional<Reindexer> node_reindexer = std::nullopt,
      const std::optional<size_t> prev_edge_count = std::nullopt,
      const std::optional<Reindexer> edge_reindexer = std::nullopt) const;

  // Build map from new NNIs to best pre-existing neighbor NNI in DAG.
  NNINNIMap BuildMapOfProposedNNIsToBestPreNNIs(const NNISet &post_nnis) const;
  // Build map from new NNI's pcsp bitsets to best reference pre-NNI's edge in DAG.
  // Creates map entry for grandparent, sister, left_child, and right_child of each NNI,
  BitsetEdgeIdMap BuildMapOfProposedNNIPCSPsToBestPreNNIEdges(
      const NNISet &post_nnis,
      std::optional<const size_t> prev_edge_count = std::nullopt,
      std::optional<const Reindexer> edge_reindexer = std::nullopt) const;
  // Build map above, but storing edges as PCSPs.
  BitsetBitsetMap BuildMapOfProposedNNIPCSPsToBestPreNNIPCSPs(
      const NNISet &post_nnis,
      std::optional<const size_t> prev_edge_count = std::nullopt,
      std::optional<const Reindexer> edge_reindexer = std::nullopt) const;

  // Build map from post-NNI PCSP to pre-NNI edge_id.
  NNIAdjBitsetEdgeIdMap BuildAdjacentPCSPsFromPreNNIToPostNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni) const;
  // Build node_id map from pre-NNI to post-NNI.
  TPChoiceMap::EdgeChoiceNodeIdMap BuildAdjacentNodeIdMapFromPreNNIToPostNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni) const;
  // Build PCSP map from pre-NNI to post-NNI.
  TPChoiceMap::EdgeChoicePCSPMap BuildAdjacentPCSPMapFromPreNNIToPostNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni) const;

  // Build map from edge PCSPs to vector of edge choice PCSPs.
  using PCSPToPCSPsMap = std::map<Bitset, std::vector<Bitset>>;
  PCSPToPCSPsMap BuildMapFromPCSPToEdgeChoicePCSPs() const;
  // Build map from edge PCSPs to their PV Hashes.
  using PCSPToPVHashesMap = std::map<Bitset, std::vector<std::string>>;
  PCSPToPVHashesMap BuildMapFromPCSPToPVHashes() const;
  // Build map from edge PCSPs to their PV Values.
  using PCSPToPVValuesMap = std::map<Bitset, std::vector<DoubleVector>>;
  PCSPToPVValuesMap BuildMapFromPCSPToPVValues() const;
  // Build map from edge PCSPs to their branch length.
  using PCSPToBranchLengthMap = std::map<Bitset, double>;
  PCSPToBranchLengthMap BuildMapFromPCSPToBranchLength() const;
  // Build map from edge PCSPs to their top tree score.
  using PCSPToScoreMap = std::map<Bitset, double>;
  PCSPToScoreMap BuildMapFromPCSPToScore(const bool recompute_scores);

  // ** TP Evaluation Engine

  // Get current in use evaluation engine.
  TPEvalEngine &GetEvalEngine() { return *eval_engine_; }
  const TPEvalEngine &GetEvalEngine() const { return *eval_engine_; }
  // Set evaluation engine type for use in runner.
  void SelectEvalEngine(const TPEvalEngineType eval_engine_type);
  // Remove all evaluation engines from use.
  void ClearEvalEngineInUse();
  // Check if evaluation engine is currently in use.
  bool IsEvalEngineInUse(const TPEvalEngineType eval_engine_type) const {
    return eval_engine_in_use_[eval_engine_type];
  }
  // Update PVs after modifying the DAG.
  void UpdateEvalEngineAfterModifyingDAG(
      const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
      const size_t prev_node_count, const Reindexer &node_reindexer,
      const size_t prev_edge_count, const Reindexer &edge_reindexer);

  // Likelihood evaluation engine.
  void MakeLikelihoodEvalEngine(const std::string &mmap_likelihood_path);
  TPEvalEngineViaLikelihood &GetLikelihoodEvalEngine() { return *likelihood_engine_; }
  const TPEvalEngineViaLikelihood &GetLikelihoodEvalEngine() const {
    return *likelihood_engine_;
  }
  bool HasLikelihoodEvalEngine() const { return likelihood_engine_ != nullptr; }
  void SelectLikelihoodEvalEngine();
  // Parsimony evaluation engine.
  void MakeParsimonyEvalEngine(const std::string &mmap_parsimony_path);
  TPEvalEngineViaParsimony &GetParsimonyEvalEngine() { return *parsimony_engine_; }
  const TPEvalEngineViaParsimony &GetParsimonyEvalEngine() const {
    return *parsimony_engine_;
  }
  bool HasParsimonyEvalEngine() const { return parsimony_engine_ != nullptr; }
  void SelectParsimonyEvalEngine();

  // Evaluation Engine - Access
  EigenConstMatrixXdRef GetLikelihoodMatrix() {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    auto &log_likelihoods =
        GetLikelihoodEvalEngine().GetDAGBranchHandler().GetBranchLengthData();
    return log_likelihoods.block(0, 0, GetEdgeCount(), log_likelihoods.cols());
  }
  PLVEdgeHandler &GetLikelihoodPVs() {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    return GetLikelihoodEvalEngine().GetPVs();
  }
  const PLVEdgeHandler &GetLikelihoodPVs() const {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    return GetLikelihoodEvalEngine().GetPVs();
  }
  PSVEdgeHandler &GetParsimonyPVs() {
    Assert(HasParsimonyEvalEngine(), "Must MakeParsimonyEvalEngine before access.");
    return GetParsimonyEvalEngine().GetPVs();
  }
  const PSVEdgeHandler &GetParsimonyPVs() const {
    Assert(HasParsimonyEvalEngine(), "Must MakeParsimonyEvalEngine before access.");
    return GetParsimonyEvalEngine().GetPVs();
  }
  EigenVectorXd &GetBranchLengths() {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    return GetLikelihoodEvalEngine().GetDAGBranchHandler().GetBranchLengthData();
  }
  const EigenVectorXd &GetBranchLengths() const {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    return GetLikelihoodEvalEngine().GetDAGBranchHandler().GetBranchLengthData();
  }
  DAGBranchHandler &GetDAGBranchHandler() {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    return GetLikelihoodEvalEngine().GetDAGBranchHandler();
  }
  const DAGBranchHandler &GetDAGBranchHandler() const {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    return GetLikelihoodEvalEngine().GetDAGBranchHandler();
  }
  EigenVectorXd &GetTopTreeLikelihoods() {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    return GetLikelihoodEvalEngine().GetTopTreeScores();
  }
  const EigenVectorXd &GetTopTreeLikelihoods() const {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    return GetLikelihoodEvalEngine().GetTopTreeScores();
  }
  EigenVectorXd &GetTopTreeParsimonies() {
    Assert(HasParsimonyEvalEngine(), "Must MakeParsimonyEvalEngine before access.");
    return GetParsimonyEvalEngine().GetTopTreeScores();
  }
  const EigenVectorXd &GetTopTreeParsimonies() const {
    Assert(HasParsimonyEvalEngine(), "Must MakeParsimonyEvalEngine before access.");
    return GetParsimonyEvalEngine().GetTopTreeScores();
  }

  // ** Branch Lengths

  // Set branch lengths from vector.
  void SetBranchLengths(EigenVectorXd new_branch_lengths) {
    GetDAGBranchHandler().SetBranchLengths(new_branch_lengths);
  }
  // Set branch lengths to default.
  void SetBranchLengthsToDefault();
  // Set branch lengths by taking the first occurrance of each PCSP edge from
  // tree collection (requires likelihood evaluation engine).
  void SetBranchLengthsByTakingFirst(const RootedTreeCollection &tree_collection,
                                     const BitsetSizeMap &edge_indexer,
                                     const bool set_uninitialized_to_default = false);
  // Find optimized branch lengths.
  void OptimizeBranchLengths(
      std::optional<bool> check_branch_convergence = std::nullopt);

  // ** Scoring

  // Get score of top-scoring tree in DAG containing given edge.
  double GetTopTreeScore(const EdgeId edge_id) const {
    return GetEvalEngine().GetTopTreeScoreWithEdge(edge_id);
  }
  // Get likelihood of top-scoring tree in DAG containing given edge.
  double GetTopTreeLikelihood(const EdgeId edge_id) const {
    Assert(HasLikelihoodEvalEngine(), "Must MakeLikelihoodEvalEngine before access.");
    return GetLikelihoodEvalEngine().GetTopTreeScoreWithEdge(edge_id);
  }
  // Get parsimony of top-scoring tree in DAG containing given edge.
  double GetTopTreeParsimony(const EdgeId edge_id) const {
    Assert(HasParsimonyEvalEngine(), "Must MakeParsimonyEvalEngine before access.");
    return GetParsimonyEvalEngine().GetTopTreeScoreWithEdge(edge_id);
  }
  // Get the Top Tree from the DAG containing the proposed NNI.
  double GetTopTreeScoreWithProposedNNI(
      const NNIOperation &post_nni, const NNIOperation &pre_nni,
      const size_t spare_offset = 0,
      std::optional<BitsetEdgeIdMap> best_edge_map_opt = std::nullopt);

  // Initialize EvalEngine.
  void InitializeScores();
  // Final Score Computation after Initialization or Update.
  void ComputeScores();
  // Update the EvalEngine after adding Node Pairs to the DAG.
  void UpdateScoresAfterDAGAddNodePair(const NNIOperation &post_nni,
                                       const NNIOperation &pre_nni,
                                       std::optional<size_t> new_tree_id);

  // ** Tree/Topology Builder

  // Get top-scoring topology in DAG containing given edge.
  Node::Topology GetTopTopologyWithEdge(const EdgeId edge_id) const;
  // Get top-scoring tree in DAG containing given edge.
  RootedTree GetTopTreeWithEdge(const EdgeId edge_id) const;

  // Get resulting top-scoring tree containing proposed NNI.
  RootedTree GetTopTreeProposedWithNNI(const NNIOperation &nni) const;
  // Get resulting top-scoring topology with proposed NNI.
  Node::Topology GetTopTreeTopologyProposedWithNNI(const NNIOperation &nni) const;

  // Build the set of edge_ids in DAG that represent the embedded top tree.
  std::set<EdgeId> BuildSetOfEdgesRepresentingTopology(
      const Node::Topology &topology) const;
  // Find the top tree's TreeId using the given edge id representation of the tree in
  // the DAG.
  std::set<TreeId> FindTreeIdsInTreeEdgeVector(const std::set<EdgeId> edge_ids) const;
  // Use branch lengths to build tree from a topology that is contained in the DAG.
  RootedTree BuildTreeFromTopologyInDAG(const Node::Topology &topology) const;

  using EdgeIdTopologyMap = std::vector<std::pair<std::set<EdgeId>, Node::Topology>>;
  using TreeIdTopologyMap = std::map<TreeId, std::vector<Node::Topology>>;
  using TreeIdTreeMap = std::map<TreeId, std::vector<RootedTree>>;
  // Build map containing all unique top tree topologies. Matched against all an
  // edge_id which results in given top tree.
  EdgeIdTopologyMap BuildMapOfEdgeIdToTopTopologies() const;
  // Build map containing all unique top tree topologies. Matched against tree_id,
  // which ranks trees according to input ordering into the DAG.
  TreeIdTopologyMap BuildMapOfTreeIdToTopTopologies() const;
  // Build map containing all unique top trees. Matched against tree_id, which ranks
  // trees according to input ordering into the DAG.
  TreeIdTreeMap BuildMapOfTreeIdToTopTrees() const;

  // Output TPEngine DAG as a newick of top trees, ordered by priority.
  std::string ToNewickOfTopTopologies() const;
  std::string ToNewickOfTopTrees() const;

  // Build PCSPs for all edges adjacent to proposed NNI.
  TPChoiceMap::EdgeChoicePCSPs BuildAdjacentPCSPsToProposedNNI(
      const NNIOperation &nni,
      const TPChoiceMap::EdgeChoiceNodeIds &adj_node_ids) const;

  // ** I/O

  std::string LikelihoodPVToString(const PVId pv_id) const;
  std::string LogLikelihoodMatrixToString() const;
  std::string ParsimonyPVToString(const PVId pv_id) const;
  std::string ChoiceMapToString() const { return GetChoiceMap().ToString(); };
  std::string TreeSourceToString() const;

 protected:
  // ** Choice Map Helpers

  // Find the edge from the highest priority tree that is adjacent to given node in
  // the given direction.
  // Accomplished by iterating over all adjacent edges using tree_source_ edge
  // map, which gives the best tree id using a given edge. The best tree is
  // expected to be the earliest found in the tree collection, aka smallest
  // tree id. The adjacent edge that comes from the best tree is chosen.
  EdgeId FindHighestPriorityEdgeAdjacentToNode(const NodeId node_id,
                                               const Direction direction) const;
  EdgeId FindHighestPriorityEdgeAdjacentToNode(const NodeId node_id,
                                               const Direction direction,
                                               const SubsplitClade clade) const;

 protected:
  // ChoiceMap for find top-scoring tree containing any given branch.
  TPChoiceMap choice_map_;
  // Tree id where branch_length and choice_map is sourced. Also function as a edge
  // priority in terms of score: lower tree_id trees should have better scores.
  std::vector<TreeId> tree_source_;
  // Total number of trees used to construct the DAG.
  size_t input_tree_count_ = 0;
  // The number of top trees in DAG.
  size_t tree_counter_ = 0;

  // Map of tree ids to topologies. The tree id gives the ranking of the best
  // scoring of inserted trees into the DAG.
  std::map<TreeId, Node::Topology> tree_id_map_;
  std::map<TreeId, double> tree_score_map_;

  // Leaf labels.
  SitePattern site_pattern_;
  EigenVectorXd site_pattern_weights_;

  size_t spare_nodes_per_nni_ = 15;
  size_t spare_edges_per_nni_ = 6;
  // Total number of nodes in DAG. Determines sizes of data vectors indexed on
  // nodes.
  size_t node_count_ = 0;
  size_t node_alloc_ = 0;
  size_t node_spare_count_ = spare_nodes_per_nni_;
  // Total number of edges in DAG. Determines sizes of data vectors indexed on edges.
  size_t edge_count_ = 0;
  size_t edge_alloc_ = 0;
  size_t edge_spare_count_ = spare_nodes_per_nni_;
  // Growth factor when reallocating data.
  constexpr static double resizing_factor_ = 2.0;

  // Un-owned reference to DAG.
  GPDAG *dag_ = nullptr;
  // Un-owned reference to GraftDAG.
  GraftDAG *graft_dag_ = nullptr;

  // Map that tracks the optimal edge to reference for each individual edge.
  BitsetEdgeIdMap best_edge_map_;
  // Use best edge map.
  bool do_use_best_edge_map_ = true;

  // A map showing which Evaluation Engines are "in use".  Several engines may be
  // instatiated, but may or may not be currently used for computation, and therefore
  // may not need to be upkept.
  TPEvalEngineTypeEnum::Array<bool> eval_engine_in_use_;
  // Un-owned reference to TP Evaluation Engine. Can be used to evaluate Top Trees
  // according to Likelihood, Parsimony, etc.
  TPEvalEngine *eval_engine_ = nullptr;
  // Engine evaluates top trees using likelihoods.
  std::unique_ptr<TPEvalEngineViaLikelihood> likelihood_engine_ = nullptr;
  // Engine evaluates top trees using parsimony.
  std::unique_ptr<TPEvalEngineViaParsimony> parsimony_engine_ = nullptr;
};
