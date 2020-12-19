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

  const SubsplitDAG dag(tree_collection_);
  auto subsplit_probabilities = DefaultDict<Bitset, double>(0.);
  const auto normalized_sbn_parameters = NormalizedSBNParameters();

  dag.IterateOverRootsplitIds(
      [&dag, &subsplit_probabilities, &normalized_sbn_parameters](size_t rootsplit_id) {
        const auto rootsplit_bitset = dag.GetDagNode(rootsplit_id)->GetBitset();
        std::cout << rootsplit_bitset.ToString() << std::endl;
        Assert(!subsplit_probabilities.contains(rootsplit_bitset),
               "We have iterated over the same rootsplit multiple times.");
        subsplit_probabilities.increment(
            rootsplit_bitset,
            normalized_sbn_parameters[dag.GetRootsplitIndex(rootsplit_bitset)]);
      });

  for (const auto node_id : dag.ReversePostorderTraversal()) {
    // TODO(e) Perhaps IterateOverEdgesAndChildren?
    const auto node = dag.GetDagNode(node_id);
    dag.IterateOverLeafwardEdgesAndNodes(
        node,
        [&node, &subsplit_probabilities, &normalized_sbn_parameters](
            const size_t gpcsp_index, const SubsplitDAGNode *child) {
          std::cout << node->GetBitset().SubsplitToString() << " ";
          std::cout << child->GetBitset().SubsplitToString() << std::endl;
          subsplit_probabilities.increment(
              child->GetBitset(), subsplit_probabilities.Map().at(node->GetBitset()) *
                                      normalized_sbn_parameters[gpcsp_index]);
        });
  }

  // Now we undo the fact that the rootsplits have been expanded in the DAG.
  // Temporary until #273.
  BitsetDoubleMap unexpanded_subsplit_probabilities;

  for (const auto &[subsplit, probability] : subsplit_probabilities.Map()) {
    auto total = subsplit.SubsplitChunk(0) | subsplit.SubsplitChunk(1);
    if (total.All()) {
      SafeInsert(unexpanded_subsplit_probabilities, subsplit.SubsplitChunk(0),
                 probability);
    } else {
      SafeInsert(unexpanded_subsplit_probabilities, subsplit, probability);
    }
  }

  return unexpanded_subsplit_probabilities;
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
