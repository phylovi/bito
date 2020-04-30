// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_instance.hpp"
#include "driver.hpp"

void GPInstance::ReadFastaFile(std::string fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

void GPInstance::ReadNewickFile(std::string fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
}

void GPInstance::ReadNexusFile(std::string fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
}

void GPInstance::CheckSequencesAndTreesLoaded() const {
  if (alignment_.SequenceCount() == 0) {
    Failwith(
        "Load an alignment into your GPInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
  if (tree_collection_.TreeCount() == 0) {
    Failwith(
        "Load some trees into your GPInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
}

void GPInstance::MakeEngine() {
  CheckSequencesAndTreesLoaded();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
  engine_ = std::make_unique<GPEngine>(site_pattern);
}

GPEngine *GPInstance::GetEngine() const {
  if (engine_ != nullptr) {
    return engine_.get();
  }
  // else
  Failwith(
      "Engine not available. Call MakeEngine to make an engine for phylogenetic "
      "likelihood computation computation.");
}

void GPInstance::ClearTreeCollectionAssociatedState() {
  sbn_parameters_.resize(0);
  rootsplits_.clear();
  indexer_.clear();
  index_to_child_.clear();
  parent_to_range_.clear();
}

void GPInstance::ProcessLoadedTrees() {
  size_t index = 0;
  ClearTreeCollectionAssociatedState();
  auto topology_counter = tree_collection_.TopologyCounter();
  // Start by adding the rootsplits.
  for (const auto &iter : RootedSBNMaps::RootsplitCounterOf(topology_counter)) {
    SafeInsert(indexer_, iter.first, index);
    rootsplits_.push_back(iter.first);
    index++;
  }
  // Now add the PCSSs.
  for (const auto &[parent, child_counter] :
       RootedSBNMaps::PCSSCounterOf(topology_counter)) {
    SafeInsert(parent_to_range_, parent, {index, index + child_counter.size()});
    for (const auto &child_iter : child_counter) {
      const auto &child = child_iter.first;
      SafeInsert(indexer_, parent + child, index);
      SafeInsert(index_to_child_, index, Bitset::ChildSubsplit(parent, child));
      index++;
    }
  }
  sbn_parameters_.resize(index);
  sbn_parameters_.setOnes();
}
