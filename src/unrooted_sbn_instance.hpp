// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "generic_sbn_instance.hpp"
#include "unrooted_sbn_support.hpp"

using PreUnrootedSBNInstance =
    GenericSBNInstance<UnrootedTreeCollection, UnrootedSBNSupport,
                       UnrootedIndexerRepresentation>;
template class GenericSBNInstance<UnrootedTreeCollection, UnrootedSBNSupport,
                                  UnrootedIndexerRepresentation>;

class UnrootedSBNInstance : public PreUnrootedSBNInstance {
 public:
  using PreUnrootedSBNInstance::PreUnrootedSBNInstance;

  // ** SBN-related items

  // max_iter is the maximum number of EM iterations to do, while score_epsilon
  // is the cutoff for score improvement.
  EigenVectorXd TrainExpectationMaximization(double alpha, size_t max_iter,
                                             double score_epsilon = 0.);

  // Sample a topology from the SBN.
  using PreUnrootedSBNInstance::SampleTopology;
  Node::NodePtr SampleTopology() const;

  // Sample trees and store them internally
  void SampleTrees(size_t count);

  // Get PSP indexer representations of the trees in tree_collection_.
  std::vector<SizeVectorVector> MakePSPIndexerRepresentations() const;

  // Return a ragged vector of vectors such that the ith vector is the
  // collection of branch lengths in the current tree collection for the ith
  // split.
  DoubleVectorVector SplitLengths() const;

  // Turn an IndexerRepresentation into a string representation of the underying
  // bitsets. This is really just so that we can make a test of indexer
  // representations.
  StringSetVector StringIndexerRepresentationOf(
      const UnrootedIndexerRepresentation &indexer_representation) const;
  StringSetVector StringIndexerRepresentationOf(const Node::NodePtr &topology,
                                                size_t out_of_sample_index) const;

  // This function is really just for testing-- it recomputes counters from
  // scratch.
  std::pair<StringSizeMap, StringPCSPMap> SplitCounters() const;

  // ** Phylogenetic likelihood

  std::vector<double> LogLikelihoods(
      std::optional<PhyloFlags> external_flags = std::nullopt);

  template <class VectorType>
  std::vector<double> LogLikelihoods(const VectorType &flag_vec,
                                     const bool is_run_defaults);

  // For each loaded tree, return the phylogenetic gradient.
  std::vector<PhyloGradient> PhyloGradients(
      std::optional<PhyloFlags> external_flags = std::nullopt);

  template <class VectorType>
  std::vector<PhyloGradient> PhyloGradients(const VectorType &flag_vec,
                                            const bool is_run_defaults);

  // Topology gradient for unrooted trees.
  // Assumption: This function is called from Python side
  // after the trees (both the topology and the branch lengths) are sampled.
  EigenVectorXd TopologyGradients(EigenVectorXdRef log_f, bool use_vimco = true);
  // Computes gradient WRT \phi of log q_{\phi}(\tau).
  // IndexerRepresentation contains all rootings of \tau.
  // normalized_sbn_parameters_in_log is a cache; see implementation of
  // TopologyGradients to see how it works. It would be private except that
  // we want to be able to test it.
  EigenVectorXd GradientOfLogQ(
      EigenVectorXdRef normalized_sbn_parameters_in_log,
      const UnrootedIndexerRepresentation &indexer_representation);

  // ** I/O

  void ReadNewickFile(const std::string &fname);
  void ReadNexusFile(const std::string &fname);

 protected:
  void PushBackRangeForParentIfAvailable(
      const Bitset &parent, UnrootedSBNInstance::RangeVector &range_vector);
  RangeVector GetSubsplitRanges(
      const RootedIndexerRepresentation &rooted_representation);
};
