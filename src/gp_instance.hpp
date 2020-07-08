// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_GP_INSTANCE_HPP_
#define SRC_GP_INSTANCE_HPP_

#include <memory.h>

#include "gp_dag.hpp"
#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "site_pattern.hpp"
#include "sugar.hpp"


size_t GetPLVIndex(PLVType plv_type, size_t node_count, size_t src_idx);

class GPInstance {
 public:
  GPInstance(std::string mmap_file_path) : mmap_file_path_(mmap_file_path) {
    if (mmap_file_path.empty()) {
      Failwith("GPInstance needs a legal path as a constructor argument.");
    }
  };
  void PrintStatus();

  void ReadFastaFile(std::string fname);
  void ReadNewickFile(std::string fname);
  void ReadNexusFile(std::string fname);

  // TODO Rename to avoid confusion with MakeEngine; we also build the DAG, so the name
  // should be more general.
  void MakeEngine();
  GPEngine *GetEngine() const;
  void PrintDAG();
  void PrintPCSPIndexer();
  void ProcessOperations(const GPOperationVector& operations);
  void EstimateSBNParameters(double tol, size_t max_iter);
  void EstimateBranchLengths(double tol, size_t max_iter);
  void PopulatePLVs();
  void ComputeLikelihoods();

  // TODO Figure out a way to remove: just here for unit test.
  size_t GetPCSPIndex(size_t parent_node_idx, size_t child_node_idx, bool rotated) {
    return dag_.GetPCSPIndex(parent_node_idx, child_node_idx, rotated);
  }

 private:
  std::string mmap_file_path_;
  Alignment alignment_;
  std::unique_ptr<GPEngine> engine_;
  RootedTreeCollection tree_collection_;
  GPDAG dag_;

  // A vector that contains all of the SBN-related probabilities.
  EigenVectorXd sbn_parameters_;
  // The master indexer for SBN parameters.

  void ClearTreeCollectionAssociatedState();
  void CheckSequencesAndTreesLoaded() const;
  void ProcessLoadedTrees();

  void MakeGPEngine();
  void InitializeGPEngine();
};

#endif  // SRC_GP_INSTANCE_HPP_
