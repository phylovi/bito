// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "gp_dag.hpp"
#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "site_pattern.hpp"
#include "nni_engine.hpp"

#include "fat_beagle.hpp"
#include "phylo_model.hpp"

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

  // ** I/O

  void ReadFastaFile(const std::string &fname);
  void ReadNewickFile(const std::string &fname, const bool sort_taxa = true);
  void ReadNewickFileGZ(const std::string &fname, const bool sort_taxa = true);
  void ReadNexusFile(const std::string &fname, const bool sort_taxa = true);
  void ReadNexusFileGZ(const std::string &fname, const bool sort_taxa = true);

  std::string GetFastaSourcePath() const;
  std::string GetNewickSourcePath() const;
  std::string GetNexusSourcePath() const;
  std::string GetMMapFilePath() const;

  void PrintStatus();
  StringSizeMap DAGSummaryStatistics();

  // ** Access

  const TagStringMap &GetTaxonMap() const;

  GPDAG &GetDAG();
  const GPDAG &GetDAG() const;

  GPEngine &GetGPEngine();
  const GPEngine &GetGPEngine() const;

  TPEngine &GetTPEngine();
  const TPEngine &GetTPEngine() const;

  NNIEngine &GetNNIEngine();
  const NNIEngine &GetNNIEngine() const;

  // ** Taxon Map

  // Get taxon names.
  StringVector GetTaxonNames() const;

  // ** DAG

  void MakeDAG();
  bool HasDAG() const;
  void PrintDAG();

  SitePattern MakeSitePattern() const;

  // Export the subsplit DAG as a DOT file.
  void SubsplitDAGToDot(const std::string &out_path,
                        bool show_index_labels = true) const;

  // ** GP Engine

  void MakeGPEngine(double rescaling_threshold = GPEngine::default_rescaling_threshold_,
                    bool use_gradients = false);
  bool HasGPEngine() const;

  void PopulatePLVs();
  void ComputeLikelihoods();
  void ComputeMarginalLikelihood();
  void CalculateHybridMarginals();

  void ResizeEngineForDAG();
  void ReinitializePriors();
  void ProcessOperations(const GPOperationVector &operations);

  void HotStartBranchLengths();
  SizeDoubleVectorMap GatherBranchLengths();
  void TakeFirstBranchLength();

  void EstimateSBNParameters();
  void SetOptimizationMethod(const OptimizationMethod method);
  void UseGradientOptimization(const bool use_gradients);
  // Estimate branch lengths using GPEngine. For testing purposes.
  void EstimateBranchLengths(double tol, size_t max_iter, bool quiet = false,
                             bool track_intermediate_iterations = false,
                             std::optional<OptimizationMethod> method = std::nullopt);

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
  RootedTreeCollection GenerateCoveringRootedTreeCollection();

  void PrintEdgeIndexer();
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

  // Get branch lengths.
  EigenVectorXd GetBranchLengths() const;
  // Get per PCSP log likelihoods
  EigenVectorXd GetPerPCSPLogLikelihoods() const;

  // ** Trees

  // Get a reference to collection of currently loaded trees.
  const RootedTreeCollection &GetCurrentlyLoadedTrees() const {
    return tree_collection_;
  };
  // Generate a version of the topologies in the current tree collection that use
  // the current GP branch lengths.
  RootedTreeCollection CurrentlyLoadedTreesWithGPBranchLengths();
  // Subset the currently loaded topologies to those that have a given PCSP, and equip
  // them with current GP branch lengths.
  RootedTreeCollection CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths(
      const std::string &pcsp_string);

  // Generate all trees spanned by the DAG and load them into the instance.
  void LoadAllGeneratedTrees();

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

  // Export spanning set of trees with GP branch lengths.
  void ExportCoveringTreesWithGPBranchLengths(const std::string &out_path) const;
  // Export spanning set of trees with TP branch lengths.
  void ExportTopTreesWithTPBranchLengths(const std::string &out_path) const;

  // ** TP Engine

  void MakeTPEngine();
  void TPEngineSetChoiceMapByTakingFirst(const bool use_subsplit_method = true);
  void TPEngineSetBranchLengthsByTakingFirst();

  // Estimate branch lengths using TPEngine. For testing purposes.
  void TPEngineEstimateBranchLengths(
      double tol, size_t max_iter, bool quiet = false,
      bool track_intermediate_iterations = false,
      std::optional<OptimizationMethod> method = std::nullopt);

  std::vector<RootedTree> TPEngineGenerateCoveringTrees();
  TreeIdTreeMap TPEngineGenerateTopRootedTrees();

  void TPEngineExportCoveringTrees(const std::string &out_path);
  void TPEngineExportTopTrees(const std::string &out_path);

  // ** NNI Evaluation Engine

  void MakeNNIEngine();

  // ** Tree Engines

  void MakeLikelihoodTreeEngine();
  FatBeagle &GetLikelihoodTreeEngine();

  void MakeParsimonyTreeEngine();
  SankoffHandler &GetParsimonyTreeEngine();

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

  // ** Data

  std::optional<std::string> fasta_path_ = std::nullopt;
  std::optional<std::string> newick_path_ = std::nullopt;
  std::optional<std::string> nexus_path_ = std::nullopt;

  RootedTreeCollection tree_collection_;
  Alignment alignment_;
  std::unique_ptr<GPDAG> dag_ = nullptr;

  // Root filepath for storing mmapped data.
  std::optional<std::string> mmap_file_path_ = std::nullopt;

  // ** Engines

  std::unique_ptr<GPEngine> gp_engine_ = nullptr;
  std::unique_ptr<TPEngine> tp_engine_ = nullptr;
  std::unique_ptr<NNIEngine> nni_engine_ = nullptr;

  std::unique_ptr<FatBeagle> likelihood_tree_engine_ = nullptr;
  std::unique_ptr<SankoffHandler> parsimony_tree_engine_ = nullptr;

  // ** Branch Length Optimization

  size_t gpcsp_count_ = 0;
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
};
