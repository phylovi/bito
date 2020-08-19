// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_sbn_instance.hpp"

StringSet RootedSBNInstance::StringIndexerRepresentationOf(
    const RootedIndexerRepresentation &indexer_representation) const {
  return RootedSBNMaps::StringIndexerRepresentationOf(
      sbn_support_.StringReversedIndexer(), indexer_representation);
}

StringSet RootedSBNInstance::StringIndexerRepresentationOf(
    const Node::NodePtr &topology, size_t out_of_sample_index) const {
  return StringIndexerRepresentationOf(
      sbn_support_.IndexerRepresentationOf(topology, out_of_sample_index));
}

std::vector<double> RootedSBNInstance::LogLikelihoods() {
  return GetEngine()->LogLikelihoods(tree_collection_, phylo_model_params_, rescaling_);
}

std::vector<RootedPhyloGradient> RootedSBNInstance::PhyloGradients() {
  return GetEngine()->Gradients(tree_collection_, phylo_model_params_, rescaling_);
}

void RootedSBNInstance::ReadNewickFile(std::string fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
  tree_collection_.ParseDatesFromTaxonNames();
  tree_collection_.InitializeParameters();
}

void RootedSBNInstance::ReadNexusFile(std::string fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
  tree_collection_.ParseDatesFromTaxonNames();
  tree_collection_.InitializeParameters();
}
