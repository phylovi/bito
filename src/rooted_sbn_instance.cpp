// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

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
  EigenVectorXd sbn_parameters = NormalizedSBNParameters();
  // Expand sbn_parameters to include fake subsplits.
  Assert(size_t(sbn_parameters.size()) == dag.EdgeCount(), "GPCSP count mismatch.");
  sbn_parameters.conservativeResize(dag.EdgeCountWithLeafSubsplits());
  sbn_parameters
      .segment(dag.EdgeCount(), dag.EdgeCountWithLeafSubsplits() - dag.EdgeCount())
      .setOnes();
  return dag.UnconditionalSubsplitProbabilities(sbn_parameters);
}

void RootedSBNInstance::UnconditionalSubsplitProbabilitiesToCSV(
    const std::string &csv_path) const {
  CSV::StringDoubleVectorToCSV(
      SBNMaps::StringDoubleVectorOf(UnconditionalSubsplitProbabilities()), csv_path);
}

std::vector<double> RootedSBNInstance::LogLikelihoods(
    std::optional<PhyloFlags> external_flags) {
  auto flags = CollectPhyloFlags(external_flags);
  return GetEngine()->LogLikelihoods(tree_collection_, phylo_model_params_, rescaling_,
                                     flags);
}

template <class VectorType>
std::vector<double> RootedSBNInstance::LogLikelihoods(const VectorType &flag_vec,
                                                      const bool is_run_defaults) {
  PhyloFlags external_flags = PhyloFlags(flag_vec, is_run_defaults);
  return LogLikelihoods(external_flags);
};
// Explicit instantiation for Pybito.
template DoubleVector RootedSBNInstance::LogLikelihoods(const StringVector &,
                                                        const bool);
template DoubleVector RootedSBNInstance::LogLikelihoods(const StringBoolVector &,
                                                        const bool);
template DoubleVector RootedSBNInstance::LogLikelihoods(const StringDoubleVector &,
                                                        const bool);
template DoubleVector RootedSBNInstance::LogLikelihoods(const StringBoolDoubleVector &,
                                                        const bool);

std::vector<double> RootedSBNInstance::UnrootedLogLikelihoods() {
  return GetEngine()->UnrootedLogLikelihoods(tree_collection_, phylo_model_params_,
                                             rescaling_);
}

std::vector<double> RootedSBNInstance::LogDetJacobianHeightTransform() {
  return GetEngine()->LogDetJacobianHeightTransform(tree_collection_,
                                                    phylo_model_params_, rescaling_);
}

std::vector<PhyloGradient> RootedSBNInstance::PhyloGradients(
    std::optional<PhyloFlags> external_flags) {
  auto flags = CollectPhyloFlags(external_flags);
  return GetEngine()->Gradients(tree_collection_, phylo_model_params_, rescaling_,
                                flags);
}

template <class VectorType>
std::vector<PhyloGradient> RootedSBNInstance::PhyloGradients(
    const VectorType &flag_vec, const bool is_run_defaults) {
  PhyloFlags external_flags = PhyloFlags(flag_vec, is_run_defaults);
  return PhyloGradients(external_flags);
};
// Explicit templates for Pybind API.
template std::vector<PhyloGradient> RootedSBNInstance::PhyloGradients(
    const StringVector &, const bool);
template std::vector<PhyloGradient> RootedSBNInstance::PhyloGradients(
    const StringBoolVector &, const bool);
template std::vector<PhyloGradient> RootedSBNInstance::PhyloGradients(
    const StringDoubleVector &, const bool);
template std::vector<PhyloGradient> RootedSBNInstance::PhyloGradients(
    const StringBoolDoubleVector &, const bool);

std::vector<DoubleVector> RootedSBNInstance::GradientLogDeterminantJacobian() {
  return GetEngine()->GradientLogDeterminantJacobian(tree_collection_,
                                                     phylo_model_params_, rescaling_);
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
