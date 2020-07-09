// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_GP_DAG_HPP_
#define SRC_GP_DAG_HPP_

#include <deque>
#include "gp_dag_node.hpp"
#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"

enum class PLVType { P, P_HAT, P_HAT_TILDE, R_HAT, R, R_TILDE };

size_t GetPLVIndex(PLVType plv_type, size_t node_count, size_t src_idx);

class GPDAG {
 public:
  GPDAG();
  explicit GPDAG(const RootedTreeCollection &tree_collection);
  void PrintStatus();

  size_t GetPCSPIndex(size_t parent_node_idx, size_t child_node_idx, bool rotated);

  size_t NodeCount() const;
  size_t ContinuousParameterCount() const;
  size_t GPCSPCount() const;
  void ProcessTrees(const RootedTreeCollection &tree_collection);
  void Print();
  void PrintPCSPIndexer();

  EigenVectorXd BuildUniformQ();

  void AddRootwardWeightedSumAccumulateOperations(std::shared_ptr<GPDAGNode> node,
                                                  bool rotated,
                                                  GPOperationVector &operations);
  void AddLeafwardWeightedSumAccumulateOperations(std::shared_ptr<GPDAGNode> node,
                                                  GPOperationVector &operations);
  void OptimizeSBNParameters(const Bitset &subsplit, GPOperationVector &operations);

  void InitializeGPEngine();
  [[nodiscard]] GPOperationVector ComputeLikelihoods();
  [[nodiscard]] GPOperationVector SetRhatToStationary();
  [[nodiscard]] GPOperationVector BranchLengthOptimization();
  [[nodiscard]] GPOperationVector SBNParameterOptimization();
  [[nodiscard]] GPOperationVector MarginalLikelihoodOperations();
  [[nodiscard]] GPOperationVector RootwardPass(std::vector<size_t> visit_order);
  [[nodiscard]] GPOperationVector RootwardPass();
  [[nodiscard]] GPOperationVector LeafwardPass(std::vector<size_t> visit_order);
  [[nodiscard]] GPOperationVector LeafwardPass();
  [[nodiscard]] GPOperationVector SetRootwardZero();
  [[nodiscard]] GPOperationVector SetLeafwardZero();
  [[nodiscard]] SizeVector LeafwardPassTraversal();
  [[nodiscard]] SizeVector RootwardPassTraversal();

  void ScheduleBranchLengthOptimization(size_t node_id,
                                        std::unordered_set<size_t> &visited_nodes,
                                        GPOperationVector &operations);
  void ScheduleSBNParametersOptimization(size_t node_id,
                                         std::unordered_set<size_t> &visited_nodes,
                                         GPOperationVector &operations);

 private:
  // TODO Can we start here with a description of the various indexing schemes? Could we
  // consider naming them different things?
  size_t taxon_count_;
  size_t gpcsp_count_;
  // A map that indexes these probabilities: rootsplits are at the beginning,
  // and PCSS bitsets are at the end.
  // The collection of rootsplits, with the same indexing as in the indexer_.
  BitsetVector rootsplits_;
  // A map going from the index of a PCSP to its child.
  SizeBitsetMap index_to_child_;
  // A map going from a parent subsplit to the range of indices in
  // sbn_parameters_ with its children. See the definition of Range for the indexing
  // convention.
  BitsetSizePairMap parent_to_range_;

  // The first entries are reserved for fake subsplits.
  // The last entries are reserved for rootsplits.
  std::unordered_map<Bitset, size_t> subsplit_to_index_;
  std::vector<std::shared_ptr<GPDAGNode>> dag_nodes_;

  // This indexer is an expanded version of indexer_ in indexer_ in sbn_instance.
  // This indexer can be used for q_, branch_lengths_, log_likelihoods_
  // in GPEngine.
  BitsetSizeMap pcsp_indexer_;
  // Stores range of indices for a subsplit and rotated subsplit.
  BitsetSizePairMap subsplit2range_;

  void CreateAndInsertNode(const Bitset &subsplit);
  void ConnectNodes(size_t idx, bool rotated);
  std::vector<Bitset> GetChildrenSubsplits(const Bitset &subsplit,
                                           bool include_fake_subsplits = false);
  void BuildNodesDepthFirst(const Bitset &subsplit,
                            std::unordered_set<Bitset> &visited_subsplits);
  void BuildNodes();
  void BuildEdges();
  void BuildPCSPIndexer();
};

#endif  // SRC_GP_DAG_HPP_
