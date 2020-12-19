// Copyright 2019-2020 libsbn project contributors.
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

  SubsplitDAG dag(tree_collection_);
  auto subsplit_probabilities = DefaultDict<Bitset, double>(0.);

  for (const auto node_id : dag.ReversePostorderTraversal()) {
    // Perhaps IterateOverEdgesAndChildren?
    dag.IterateOverLeafwardEdgesAndNodes(
        dag.GetDagNode(node_id),
        [this, &subsplit_probabilities](const size_t gpcsp_index,
                                        const SubsplitDAGNode *child) {
          subsplit_probabilities.increment(child->GetBitset(),
                                           sbn_parameters_[gpcsp_index]);
        });
  }

  // TODO(e) add the rootsplits.

  return subsplit_probabilities.Map();
}

std::vector<double> RootedSBNInstance::LogLikelihoods() {
  return GetEngine()->LogLikelihoods(tree_collection_, phylo_model_params_, rescaling_);
}

std::vector<RootedPhyloGradient> RootedSBNInstance::PhyloGradients() {
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
