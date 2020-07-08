// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_GP_DAG_HPP_
#define SRC_GP_DAG_HPP_

#include <memory.h>

#include <deque>

#include "dag_node.hpp"
#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "site_pattern.hpp"
#include "sugar.hpp"

enum class PLVType { P, P_HAT, P_HAT_TILDE, R_HAT, R, R_TILDE };

size_t GetPLVIndex(PLVType plv_type, size_t node_count, size_t src_idx);

class GPDAG {
 public:
  GPDAG();
  explicit GPDAG(const RootedTreeCollection &tree_collection);
  void PrintStatus();

  size_t GetPCSPIndex(size_t parent_node_idx, size_t child_node_idx, bool rotated);

  size_t NodeCount() const;
  size_t PCSPCount() const;
  size_t GPCSPCount() const;
  size_t SBNParameterCount() const;
  void ProcessTrees(const RootedTreeCollection &tree_collection);
  void Construct();
  void Print();

  EigenVectorXd BuildUniformQ();

  GPOperationVector ComputeLikelihoods();
  GPOperationVector SetRhatToStationary();
  void AddRootwardWeightedSumAccumulateOperations(std::shared_ptr<DAGNode> node,
                                                  bool rotated,
                                                  GPOperationVector &operations);
  void AddLeafwardWeightedSumAccumulateOperations(std::shared_ptr<DAGNode> node,
                                                  GPOperationVector &operations);
  void OptimizeSBNParameters(const Bitset &subsplit, GPOperationVector &operations);

  void InitializeGPEngine();
  GPOperationVector SetStationaryDistribution();
  GPOperationVector BranchLengthOptimization();
  GPOperationVector SBNParameterOptimization();
  GPOperationVector MarginalLikelihoodOperations();
  void ScheduleBranchLengthOptimization(size_t node_id,
                                        std::unordered_set<size_t> &visited_nodes,
                                        GPOperationVector &operations);
  void ScheduleSBNParametersOptimization(size_t node_id,
                                         std::unordered_set<size_t> &visited_nodes,
                                         GPOperationVector &operations);
  GPOperationVector RootwardPass(std::vector<size_t> visit_order);
  GPOperationVector RootwardPass();
  GPOperationVector LeafwardPass(std::vector<size_t> visit_order);
  GPOperationVector LeafwardPass();
  GPOperationVector SetRootwardZero();
  GPOperationVector SetLeafwardZero();
  std::vector<size_t> LeafwardPassTraversal();
  std::vector<size_t> RootwardPassTraversal();

 private:
  size_t taxon_count_;
  size_t sbn_parameter_count_;
  // The master indexer for SBN parameters.
  BitsetSizeMap indexer_;
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
  // We will adopt the convention of denoting a trivial subsplit by placing
  // 0's in the first chunk.
  // The next set of entries store the rootsplits.
  // The remaining entries - 1 are for all other subsplits.
  // The last entry is saved for pseudoroot -- for passing stationary.
  // The pseudoroot is a Bitset with all 1's.
  std::unordered_map<Bitset, size_t> subsplit_to_index_;
  std::vector<std::shared_ptr<DAGNode>> dag_nodes_;

  // Build indexing scheme for q_, branch_lengths_, log_likelihoods_.
  // PCSP to index.
  BitsetSizeMap pcsp_indexer_;
  // Stores range of indices for a subsplit and rotated subsplit.
  BitsetSizePairMap subsplit2range_;
  std::vector<size_t> rootward_order_;
  std::vector<size_t> leafward_order_;

  void CreateAndInsertNode(const Bitset &subsplit);
  void ConnectNodes(size_t idx, bool rotated);
  void AddChildrenSubsplits(const Bitset &subsplit,
                            std::deque<Bitset> &subsplit_queue,
                            std::unordered_set<Bitset> &visited_subsplits);
  std::vector<Bitset> GetChildrenSubsplits(const Bitset &subsplit,
                                           bool include_fake_subsplits = false);
  void BuildNodesDepthFirst(const Bitset &subsplit,
                            std::deque<Bitset> &subsplit_queue,
                            std::unordered_set<Bitset> &visited_subsplits);
  void BuildNodes();
  void BuildEdges();
  void BuildPCSPIndexer();
};

#endif  // SRC_GP_DAG_HPP_
