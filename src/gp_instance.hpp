// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "gp_dag.hpp"
#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "site_pattern.hpp"
#include "nni_engine.hpp"

// New typedef used for storing/outputting intermediate or perturbed+tracked values
// from branch length estimation.
using VectorOfStringAndEigenVectorXdPairs =
    std::vector<std::pair<std::string, EigenVectorXd>>;

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
  void UseGradientOptimization(bool use_gradients = false);

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
                             bool track_intermediate_iterations = false);
  void PopulatePLVs();
  void ComputeLikelihoods();
  void ComputeMarginalLikelihood();
  void CalculateHybridMarginals();

  // This scans the PCSP likelihood surface by calculating the per pcsp likelihood
  // values at different branch length values. The currently set branch lengths are
  // scaled by a vector of size "steps" that ranges linearly from "scale_min" to
  // "scale_max".
  void GetPerGPCSPLogLikelihoodSurfaces(size_t steps, double scale_min,
                                        double scale_max);
  // This is for tracking branch length optimization following perturbation of a single
  // branch length, when assuming all other branch lengths are optimal. We perturb
  // branch lengths for each pcsp to the default value of 0.1 and then track branch
  // length and per pcsp likelihood values until the likelihood converges to the optimal
  // value or the number of DAG traversals exceeds 5.
  void PerturbAndTrackValuesFromOptimization();
  RootedTreeCollection GenerateCompleteRootedTreeCollection();

  // #348: A lot of code duplication here with things in SBNInstance.
  StringVector PrettyIndexer() const;
  EigenConstVectorXdRef GetSBNParameters();
  StringDoubleVector PrettyIndexedSBNParameters();
  StringDoubleVector PrettyIndexedBranchLengths();
  StringDoubleVector PrettyIndexedPerGPCSPLogLikelihoods();
  StringDoubleVector PrettyIndexedPerGPCSPComponentsOfFullLogMarginal();

  VectorOfStringAndEigenVectorXdPairs PrettyIndexedIntermediateBranchLengths();
  VectorOfStringAndEigenVectorXdPairs PrettyIndexedIntermediatePerGPCSPLogLikelihoods();
  VectorOfStringAndEigenVectorXdPairs PrettyIndexedPerGPCSPLogLikelihoodSurfaces();

  void SBNParametersToCSV(const std::string &file_path);
  void SBNPriorToCSV(const std::string &file_path);
  void BranchLengthsToCSV(const std::string &file_path);
  void PerGPCSPLogLikelihoodsToCSV(const std::string &file_path);
  void IntermediateBranchLengthsToCSV(const std::string &file_path);
  void IntermediatePerGPCSPLogLikelihoodsToCSV(const std::string &file_path);
  void PerGPCSPLogLikelihoodSurfacesToCSV(const std::string &file_path);
  void TrackedOptimizationValuesToCSV(const std::string &file_path);

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

  // Calculate and store the intermediate per pcsp branch length and likelihood values
  // during branch length estimation, so that they can be output to CSV.
  void IntermediateOptimizationValues();

  EdgeId GetEdgeIndexForLeafNode(const Bitset &parent_subsplit,
                                 const Node *leaf_node) const;
  RootedTreeCollection TreesWithGPBranchLengthsOfTopologies(
      Node::NodePtrVec &&topologies) const;
  StringDoubleVector PrettyIndexedVector(EigenConstVectorXdRef v);
  VectorOfStringAndEigenVectorXdPairs PrettyIndexedMatrix(EigenConstMatrixXdRef m);
  void PerPCSPIndexedMatrixToCSV(
      VectorOfStringAndEigenVectorXdPairs per_pcsp_indexed_matrix,
      const std::string &file_path);

  std::string mmap_file_path_;
  Alignment alignment_;
  std::unique_ptr<GPEngine> engine_;
  RootedTreeCollection tree_collection_;
  GPDAG dag_;
  static constexpr size_t plv_count_per_node_ = 6;

  size_t gpcsp_count_ = dag_.EdgeCountWithLeafSubsplits();
  bool use_gradients_;
  GPEngine::OptimizationMethod optimization_method_;

  // For storing intermediate optimization branch length and per pcsp log
  // likelihood values. Only used if track_intermediate_iterations in
  // EstimateBranchLengths is true.
  EigenMatrixXd per_pcsp_branch_lengths_ = EigenMatrixXd(gpcsp_count_, 1);
  EigenMatrixXd per_pcsp_log_lik_ = EigenMatrixXd(gpcsp_count_, 1);
  // For storing branch length and log likelihood values when finding the log likelihood
  // surface for each pcsp
  EigenMatrixXd per_pcsp_lik_surfaces_;
  // For storing outputs after perturbing and then tracking branch length and per pcsp
  // log likelihoods
  VectorOfStringAndEigenVectorXdPairs tracked_values_after_perturbing_;

  std::unique_ptr<NNIEngine> nni_engine_;
};
