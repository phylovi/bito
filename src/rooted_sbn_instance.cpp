// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_sbn_instance.hpp"
#include "subsplit_dag.hpp"

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

BitsetDoubleMap RootedSBNInstance::UnconditionalSubsplitProbabilities() const {
  if (tree_collection_.TreeCount() == 0) {
    Failwith(
        "Please load some trees into your RootedSBNInstance before trying to calculate "
        "UnconditionalSubsplitProbabilities.");
  }
  const SubsplitDAG dag(tree_collection_);
  return dag.UnconditionalSubsplitProbabilities(NormalizedSBNParameters());
}

void RootedSBNInstance::UnconditionalSubsplitProbabilitiesToCSV(
    const std::string &csv_path) const {
  CSV::StringDoubleVectorToCSV(
      SBNMaps::StringDoubleVectorOf(UnconditionalSubsplitProbabilities()), csv_path);
}

std::vector<double> RootedSBNInstance::LogLikelihoods() {
  return GetEngine()->LogLikelihoods(tree_collection_, phylo_model_params_, rescaling_);
}

std::vector<double> RootedSBNInstance::UnrootedLogLikelihoods() {
  return GetEngine()->UnrootedLogLikelihoods(tree_collection_, phylo_model_params_,
                                             rescaling_);
}

std::vector<PhyloGradient> RootedSBNInstance::PhyloGradients() {
  return GetEngine()->Gradients(tree_collection_, phylo_model_params_, rescaling_);
}

void RootedSBNInstance::ReadNewickFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
}

void RootedSBNInstance::ReadNexusFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
}

void RootedSBNInstance::SetDatesToBeConstant(
    bool initialize_time_trees_using_branch_lengths) {
  tree_collection_.SetDatesToBeConstant(initialize_time_trees_using_branch_lengths);
}

void RootedSBNInstance::ParseDatesFromTaxonNames(
    bool initialize_time_trees_using_branch_lengths) {
  tree_collection_.ParseDatesFromTaxonNames(initialize_time_trees_using_branch_lengths);
}

void RootedSBNInstance::ParseDatesFromCSV(
    const std::string &csv_path, bool initialize_time_trees_using_branch_lengths) {
  tree_collection_.ParseDatesFromCSV(csv_path,
                                     initialize_time_trees_using_branch_lengths);
}
