// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "gp_dag.hpp"
#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "site_pattern.hpp"
#include "nni_engine.hpp"

// New typedef
using StringEigenVectorXdVector = std::vector<std::pair<std::string, EigenVectorXd>>;

class GPInstance {
 public:
  explicit GPInstance(const std::string &mmap_file_path)
      : mmap_file_path_(mmap_file_path) {
    if (mmap_file_path.empty()) {
      Failwith("GPInstance needs a legal path as a constructor argument.");
    }
  };
  void PrintStatus();
  StringSizeMap DAGSummaryStatistics();

  void ReadFastaFile(const std::string &fname);
  void ReadNewickFile(const std::string &fname);
  void ReadNewickFileGZ(const std::string &fname);
  void ReadNexusFile(const std::string &fname);
  void ReadNexusFileGZ(const std::string &fname);

  void MakeDAG();
  GPDAG &GetDAG();
  void PrintDAG();

  void MakeEngine(double rescaling_threshold = GPEngine::default_rescaling_threshold_);
  GPEngine *GetEngine() const;
  bool HasEngine() const;
  void ResizeEngineForDAG();

  void PrintEdgeIndexer();
  void ReinitializePriors();
  void ProcessOperations(const GPOperationVector &operations);
  void HotStartBranchLengths();
  SizeDoubleVectorMap GatherBranchLengths();
  void TakeFirstBranchLength();
  void EstimateSBNParameters();
  void EstimateBranchLengths(double tol, size_t max_iter, bool quiet = false,
                             int optim_tol = 10);
  void PopulatePLVs();
  void ComputeLikelihoods();
  void ComputeMarginalLikelihood();
  void CalculateHybridMarginals();
  void GetPerGPCSPLogLikelihoodSurfaces(int steps);
  RootedTreeCollection GenerateCompleteRootedTreeCollection();

  // #348: A lot of code duplication here with things in SBNInstance.
  StringVector PrettyIndexer() const;
  EigenConstVectorXdRef GetSBNParameters();
  StringDoubleVector PrettyIndexedSBNParameters();
  StringDoubleVector PrettyIndexedBranchLengths();
  StringDoubleVector PrettyIndexedPerGPCSPLogLikelihoods();
  StringDoubleVector PrettyIndexedPerGPCSPComponentsOfFullLogMarginal();

  StringEigenVectorXdVector PrettyIndexedPerGPCSPBranchLengthsFromOptimization();
  StringEigenVectorXdVector PrettyIndexedPerGPCSPLogLikelihoodsFromOptimization();
  StringEigenVectorXdVector PrettyIndexedPerGPCSPLogLikelihoodSurfaces();
  void TrackValuesFromOptimization();

  void SBNParametersToCSV(const std::string &file_path);
  void SBNPriorToCSV(const std::string &file_path);
  void BranchLengthsToCSV(const std::string &file_path);
  void PerGPCSPLogLikelihoodsToCSV(const std::string &file_path);
  void PerGPCSPBranchLengthsFromOptimizationToCSV(const std::string &file_path);
  void PerGPCSPLogLikelihoodsFromOptimizationToCSV(const std::string &file_path);
  void PerGPCSPLogLikelihoodSurfacesToCSV(const std::string &file_path);
  void FullDAGTraversalOptimizationValuesToCSV(const std::string &file_path);

  // Generate a version of the topologies in the current tree collection that use
  // the current GP branch lengths.
  RootedTreeCollection CurrentlyLoadedTreesWithGPBranchLengths();

  // Subset the currently loaded topologies to those that have a given PCSP, and equip
  // them with current GP branch lengths.
  RootedTreeCollection CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths(
      const std::string &pcsp_string);

  // Run CurrentlyLoadedTreesWithGPBranchLengths and export to a Newick file.
  void ExportTrees(const std::string &out_path);
  // Run CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths and export to a Newick
  // file.
  void ExportTreesWithAPCSP(const std::string &pcsp_string,
                            const std::string &newick_path);
  // Export all topologies in the span of the subsplit DAG to a Newick file. Does not
  // require an Engine.
  void ExportAllGeneratedTopologies(const std::string &out_path);
  // Export all trees in the span of the subsplit DAG (with GP branch lengths) to a
  // Newick file. Requires an Engine.
  void ExportAllGeneratedTrees(const std::string &out_path);
  // Generate all trees spanned by the DAG and load them into the instance.
  void LoadAllGeneratedTrees();

  // Get taxon names.
  StringVector GetTaxonNames() const;
  // Get branch lengths.
  EigenVectorXd GetBranchLengths() const;
  // Export the subsplit DAG as a DOT file.
  void SubsplitDAGToDot(const std::string &out_path,
                        bool show_index_labels = true) const;

  // Initialize NNI Evaluation Engine.
  void MakeNNIEngine();
  // Get NNI Evaluation Engine.
  NNIEngine &GetNNIEngine();

 private:
  void ClearTreeCollectionAssociatedState();
  void CheckSequencesLoaded() const;
  void CheckTreesLoaded() const;

  size_t GetEdgeIndexForLeafNode(const Bitset &parent_subsplit,
                                 const Node *leaf_node) const;
  RootedTreeCollection TreesWithGPBranchLengthsOfTopologies(
      Node::NodePtrVec &&topologies) const;

  std::string mmap_file_path_;
  Alignment alignment_;
  std::unique_ptr<GPEngine> engine_;
  RootedTreeCollection tree_collection_;
  GPDAG dag_;
  static constexpr size_t plv_count_per_node_ = 6;

  // For optimization tracking
  StringEigenVectorXdVector tracked_values_after_full_dag_traversal_;
};
