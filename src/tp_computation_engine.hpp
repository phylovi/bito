// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// TP Computation Engine

#pragma once

#include "tp_engine.hpp"

class TPComputationEngine {
 public:
  // Initialize Computation Engine.
  virtual void Init();
  // Update the Computation Engine after adding Node Pairs to the DAG.
  virtual void UpdateAfterDAGNodePair(const NNIOperation &post_nni,
                                      const NNIOperation &pre_nni,
                                      std::optional<size_t> new_tree_id);
  // Get the Top Tree from the DAG with the given edge.
  virtual double GetTopTreeLikelihoodWithEdge(const EdgeId edge_id);
};

class TPComputationEngineViaLikelihood : public TPComputationEngine {
 public:
  // Initialize Computation Engine.
  virtual void Init();
  // Update the Computation Engine after adding Node Pairs to the DAG.
  virtual void UpdateAfterDAGNodePair(const NNIOperation &post_nni,
                                      const NNIOperation &pre_nni,
                                      std::optional<size_t> new_tree_id);
  // Get the Top Tree from the DAG with the given edge.
  virtual double GetTopTreeLikelihoodWithEdge(const EdgeId edge_id);

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

  // Tree likelihoods matrix across all sites.
  EigenMatrixXd log_likelihoods_;
  // Top tree log likelihood per edge.
  EigenVectorXd top_tree_log_likelihoods_per_edge_;
};

class TPComputationEngineViaParsimony : public TPComputationEngine {
 public:
  // Initialize Computation Engine.
  virtual void Init();
  // Update the Computation Engine after adding Node Pairs to the DAG.
  virtual void UpdateAfterDAGNodePair(const NNIOperation &post_nni,
                                      const NNIOperation &pre_nni,
                                      std::optional<size_t> new_tree_id);
  // Get the Top Tree from the DAG with the given edge.
  virtual double GetTopTreeLikelihoodWithEdge(const EdgeId edge_id);

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

  // Tree parsimonies across all sites.
  EigenMatrixXd parsimonies_;
  // Top tree parsimony per edge.
  EigenVectorXd top_tree_parsimony_per_edge_;
};
