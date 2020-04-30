// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_GP_INSTANCE_HPP_
#define SRC_GP_INSTANCE_HPP_

#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "site_pattern.hpp"
#include "sugar.hpp"

class GPInstance {
 public:
  GPInstance(){};

  void ReadFastaFile(std::string fname);
  void ReadNewickFile(std::string fname);
  void ReadNexusFile(std::string fname);

  void MakeEngine();

 private:
  Alignment alignment_;
  std::unique_ptr<GPEngine> engine_;
  RootedTreeCollection tree_collection_;

  // A vector that contains all of the SBN-related probabilities.
  EigenVectorXd sbn_parameters_;
  // The master indexer for SBN parameters.
  BitsetSizeMap indexer_;
  // A vector of the taxon names.
  StringVector taxon_names_;
  // A map that indexes these probabilities: rootsplits are at the beginning,
  // and PCSS bitsets are at the end.
  // The collection of rootsplits, with the same indexing as in the indexer_.
  BitsetVector rootsplits_;
  // A map going from the index of a PCSS to its child.
  SizeBitsetMap index_to_child_;
  // A map going from a parent subsplit to the range of indices in
  // sbn_parameters_ with its children. See the definition of Range for the indexing
  // convention.
  BitsetSizePairMap parent_to_range_;

  void ClearTreeCollectionAssociatedState();
  void CheckSequencesAndTreesLoaded() const;
  void ProcessLoadedTrees();
  GPEngine *GetEngine() const;
};




#endif  // SRC_GP_INSTANCE_HPP_
