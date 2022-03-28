// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// NNI Engine
// This engine is used to explore the topological space surrounding the
// subsplitDAG using the smallest topological change of the NNI (Nearest Neighbor
// Interchange). The engine functions by finding all the NNIs adjacent to the DAG,
// evaluating them by their likelihoods, then choosing to accept or reject each NNI
// based on a filtering criteria.  For each NNI added to the DAG, each new edge added
// also adds new adjacent NNIs to the DAG.  This process is repeated until no NNIs pass
// the filtering criteria.

#pragma once

#include "gp_engine.hpp"
#include "gp_dag.hpp"

#include "bitset.hpp"
#include "subsplit_dag.hpp"
#include "nni_operation.hpp"
#include "graft_dag.hpp"
#include "sugar.hpp"

class NNIEngine {
 public:
  // Constructors:
  NNIEngine(GPDAG &dag, GPEngine &gp_engine);

  // ** Getters

  // Get Reference of DAG.
  const GPDAG &GetGPDAG() const { return dag_; };
  // Get Reference of GraftDAG.
  GraftDAG &GetGraftDAG() const { return *graft_dag_.get(); };
  // Get Reference of GPEngine.
  const GPEngine &GetGPEngine() const { return gp_engine_; }
  // Get Adjacent NNIs to DAG.
  const NNISet &GetAdjacentNNIs() const { return adjacent_nnis_; };
  // Get Number of Adjacent NNIs.
  size_t GetAdjacentNNICount() const { return adjacent_nnis_.size(); };
  // Get All NNIs that have been accepted into the DAG.
  const std::vector<NNIOperation> &GetAcceptedNNIs() const { return accepted_nnis_; };
  // Get All NNIs that have been rejected from the DAG.
  const std::vector<NNIOperation> &GetRejectedNNIs() const { return rejected_nnis_; };
  // Get Number of Accepted NNIs.
  size_t GetAcceptedNNICount() const { return accepted_nnis_.size(); };
  // Get Map of Adjacent NNIs with their score.
  const std::map<NNIOperation, double> &GetScoredNNIs() const { return scored_nnis_; };
  // Get number of runs of NNI engine.
  size_t GetRunCount() const { return run_count_; };

  // ** Runners
  // These start the engine, which procedurally ranks and adds (and
  // maybe removes) NNIs to the DAG, until some termination criteria has been satisfied.

  // Primary Runner using GraftDAG.
  void Run();
  // Runner subroutines.
  void RunInit();
  void RunMainLoop();
  void RunPostLoop();

  // ** Evaluation/Ranking
  // These evaluate all NNIs from the NNISet by a criterion in order to place
  // them in a relative ordering.

  // Evaluate all adjacent NNIs.
  void EvaluateAdjacentNNIs();
  // Method Selector for Individual NNI Evaluation.
  double EvaluateNNI(const NNIOperation &proposed_nni);
  // Compute NNI Likelihood, via GPEngine's quartet request.
  double ComputeNNILikelihood(const NNIOperation &nni);
  // Add score to given NNI.
  void AddScoreForNNI(const NNIOperation &nni, const double score);
  double GetScoreForNNI(const NNIOperation &nni) const;

  // ** Filtering

  // Function template for filtering adjacent NNI to accepted or rejected.
  using FilterFunction =
      std::function<bool(NNIEngine &, GPEngine &, const NNIOperation &)>;
  static bool NoFilter(NNIEngine &this_nni_engine, GPEngine &this_gp_engine,
                       const NNIOperation &nni);

  // Initialize filter before first NNI sweep.
  void FilterInit();
  // Update filter parameters for each NNI sweep.
  void FilterUpdate();
  // Apply the filtering method to determine whether each Adjacent NNI will be accepted
  // or rejected.
  // Note: default function accepts all.
  void FilterAdjacentNNIs();

  // ** DAG & GPEngine Maintenance

  // Initial GPEngine for use with GraftDAG.
  void InitGPEngine();
  // Resize GPengine for GraftDAG current size.
  void ResizeGPEngineForDAG();
  // Resize GPengine for GraftDAG current size.
  void ResizeGPEngineForGraftDAG();
  // Populate PLVs so that
  void PrepGPEngineForLikelihoods();

  // Add all Accepted NNIs to Main DAG.
  void AddAcceptedNNIsToDAG();
  // Add all Adjacent NNIs to Graft DAG.
  void GraftAdjacentNNIsToDAG();
  // Remove all NNIs from Graft DAG.
  void RemoveAllGraftedNNIsFromDAG();

  // ** NNI Maintenance
  // These maintain NNIs to stay consistent with the state of associated GraftDAG.

  // Freshly synchonizes NNISet to match the current state of its DAG. Wipes old NNI
  // data and finds all all parent/child pairs adjacent to DAG by iterating over all
  // internal edges in the DAG. (For each internal edges, two NNIs are possible.)
  void SyncAdjacentNNIsWithDAG();
  // Updates NNI Set after given parent/child node pair have been added to the DAG.
  // Removes pair from NNI Set and adds adjacent pairs coming from newly created edges.
  void UpdateAdjacentNNIsAfterDAGAddNodePair(const NNIOperation &nni);
  void UpdateAdjacentNNIsAfterDAGAddNodePair(const Bitset &parent_bitset,
                                             const Bitset &child_bitset);
  // Adds all adjacent NNIs from new created edges by accepted NNIs.
  void UpdateAdjacentNNIsAfterAddAcceptedNNIs();
  // Adds all NNIs from all (node_id, other_id) pairs, where other_id's are elements of
  // the adjacent_node_ids vector. is_edge_leafward tells whether node_id is the child
  // or parent. is_edge_edgeclade determines which side of parent the child descends
  // from.
  void AddAllNNIsFromNodeVectorToAdjacentNNIs(const size_t &node_id,
                                              const SizeVector &adjacent_node_ids,
                                              const bool is_edge_on_left,
                                              const bool is_edge_leafward);
  // Based on given input NNIOperation, produces the two possible output NNIOperations
  // and adds those results to the NNI Set (if results are not a member of the DAG).
  void SafeAddOutputNNIsToAdjacentNNIs(const Bitset &parent_bitset,
                                       const Bitset &child_bitset,
                                       const bool is_edge_on_left);
  // Remove all adjacent nnis from list.
  void ClearAdjacentNNIs();
  // Remove all accepted nnis from list.
  void ClearAcceptedNNIs();
  // Remove all rejected nnis from list.
  void ClearRejectedNNIs();

 private:
  GPDAG &dag_;
  std::unique_ptr<GraftDAG> graft_dag_;
  GPEngine &gp_engine_;

  // Set of NNIs to be evaluated, which are a single NNI.
  NNISet adjacent_nnis_;
  // NNIs which have passed the filtering threshold, to be added to the DAG.
  std::vector<NNIOperation> accepted_nnis_;
  // NNIs which have failed the filtering threshold, to be added to the DAG.
  std::vector<NNIOperation> rejected_nnis_;
  // Map of adjacent NNIs to their score.
  std::map<NNIOperation, double> scored_nnis_;
  // Count number of loops executed by engine.
  size_t run_count_;
};
