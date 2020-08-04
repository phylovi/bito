// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_sbn_instance.hpp"

void RootedSBNInstance::TrainSimpleAverage() {
  CheckTopologyCounter();
  auto indexer_representation_counter = RootedSBNMaps::IndexerRepresentationCounterOf(
      indexer_, topology_counter_, sbn_parameters_.size());
  SBNProbability::SimpleAverage(sbn_parameters_, indexer_representation_counter,
                                rootsplits_.size(), parent_to_range_);
}

void RootedSBNInstance::ReadNewickFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
  tree_collection_.ParseDatesFromTaxonNames();
  tree_collection_.InitializeParameters();
}

void RootedSBNInstance::ReadNexusFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
  tree_collection_.ParseDatesFromTaxonNames();
  tree_collection_.InitializeParameters();
}

std::vector<double> RootedSBNInstance::LogLikelihoods() {
  return GetEngine()->LogLikelihoods(tree_collection_, phylo_model_params_, rescaling_);
}

std::vector<RootedPhyloGradient> RootedSBNInstance::PhyloGradients() {
  return GetEngine()->Gradients(tree_collection_, phylo_model_params_, rescaling_);
}


