// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "csv.hpp"
#include "generic_sbn_instance.hpp"
#include "phylo_flags.hpp"
#include "rooted_gradient_transforms.hpp"
#include "rooted_sbn_support.hpp"

using PreRootedSBNInstance = GenericSBNInstance<RootedTreeCollection, RootedSBNSupport,
                                                RootedIndexerRepresentation>;
template class GenericSBNInstance<RootedTreeCollection, RootedSBNSupport,
                                  RootedIndexerRepresentation>;

class RootedSBNInstance : public PreRootedSBNInstance {
 public:
  using PreRootedSBNInstance::PreRootedSBNInstance;

  // ** SBN-related items

  // Turn an IndexerRepresentation into a string representation of the underlying
  // bitsets. This is really just so that we can make a test of indexer
  // representations.
  StringSet StringIndexerRepresentationOf(
      const RootedIndexerRepresentation& indexer_representation) const;
  StringSet StringIndexerRepresentationOf(const Node::NodePtr& topology,
                                          size_t out_of_sample_index) const;

  // Make a map from each subsplit to its overall probability when we sample a tree from
  // the SBN.
  BitsetDoubleMap UnconditionalSubsplitProbabilities() const;
  void UnconditionalSubsplitProbabilitiesToCSV(const std::string& csv_path) const;

  // ** Phylogenetic likelihood

  std::vector<double> LogLikelihoods(
      std::optional<PhyloFlags> external_flags = std::nullopt);

  template <class VectorType>
  std::vector<double> LogLikelihoods(const VectorType& flag_vec,
                                     const bool is_run_defaults);

  std::vector<double> UnrootedLogLikelihoods();

  std::vector<double> LogDetJacobianHeightTransform();

  std::vector<PhyloGradient> PhyloGradients(
      std::optional<PhyloFlags> external_flags = std::nullopt);

  template <class VectorType>
  std::vector<PhyloGradient> PhyloGradients(const VectorType& flag_vec,
                                            const bool is_run_defaults);

  std::vector<DoubleVector> GradientLogDeterminantJacobian();

  // ** I/O

  void ReadNewickFile(const std::string& fname);
  void ReadNexusFile(const std::string& fname);

  void SetDatesToBeConstant(bool initialize_time_trees_using_branch_lengths);
  void ParseDatesFromTaxonNames(bool initialize_time_trees_using_branch_lengths);
  void ParseDatesFromCSV(const std::string& csv_path,
                         bool initialize_time_trees_using_branch_lengths);
};
