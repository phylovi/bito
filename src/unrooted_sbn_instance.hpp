// Copyright 2019-2021 bito project contributors.
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

  std::vector<double> LogLikelihoods();

  // For each loaded tree, return the phylogenetic gradient.
  std::vector<PhyloGradient> PhyloGradients();
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

#ifdef DOCTEST_LIBRARY_INCLUDED

#include "doctest_constants.hpp"

TEST_CASE("UnrootedSBNInstance: indexer and PSP representations") {
  UnrootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon_unrooted.nwk");
  inst.ProcessLoadedTrees();
  auto pretty_indexer = inst.PrettyIndexer();
  // The indexer_ is to index the sbn_parameters_. Note that neither of these
  // data structures attempt to catalog the complete collection of rootsplits or
  // PCSPs, but just those that are present for some rooting of the input trees.
  //
  // The indexer_ and sbn_parameters_ are laid out as follows (I'll just call it
  // the "index" in what follows). Say there are rootsplit_count rootsplits in
  // the support.
  // The first rootsplit_count entries of the index are assigned to the
  // rootsplits (again, those rootsplits that are present for some rooting of
  // the unrooted input trees). For the five_taxon example, this goes as follows:
  StringSet correct_pretty_rootsplits(
      {"00000|11111|01110", "00000|11111|01010", "00000|11111|00101",
       "00000|11111|00111", "00000|11111|00001", "00000|11111|00011",
       "00000|11111|00010", "00000|11111|00100", "00000|11111|00110",
       "00000|11111|01000", "00000|11111|01111", "00000|11111|01001"});
  StringSet pretty_rootsplits(
      pretty_indexer.begin(),
      pretty_indexer.begin() + correct_pretty_rootsplits.size());
  CHECK(correct_pretty_rootsplits == pretty_rootsplits);
  // The rest of the entries of the index are laid out as blocks of parameters
  // for PCSPs that share the same parent. Take a look at the description of
  // PCSP bitsets (and the unit tests) in bitset.hpp to understand the notation
  // used here.
  //
  // For example, here are four PCSPs that all share the parent 00001|11110:
  StringSet correct_pretty_pcsp_block({"00001|11110|01110", "00001|11110|00010",
                                       "00001|11110|01000", "00001|11110|00100"});
  StringSet pretty_indexer_set(pretty_indexer.begin(), pretty_indexer.end());
  // It's true that this test doesn't show the block-ness, but it wasn't easy to
  // show off this feature in a way that wasn't compiler dependent.
  // You can see it by printing out a pretty_indexer if you wish. A test exhibiting
  // block structure appeas in rooted_sbn_instance.hpp.
  for (auto pretty_pcsp : correct_pretty_pcsp_block) {
    CHECK(pretty_indexer_set.find(pretty_pcsp) != pretty_indexer_set.end());
  }
  // Now we can look at some tree representations. We get these by calling
  // IndexerRepresentationOf on a tree topology. This function "digests" the
  // tree by representing all of the PCSPs as bitsets which it can then look up
  // in the indexer_.
  // It then spits them out as the rootsplit and PCSP indices.
  // The following tree is (2,(1,3),(0,4));, or with internal nodes (2,(1,3)5,(0,4)6)7
  auto indexer_test_topology_1 = Node::OfParentIdVector({6, 5, 7, 5, 6, 7, 7});
  // Here we look at the indexer representation of this tree. Rather than having
  // the indices themselves, which is what IndexerRepresentationOf actually
  // outputs, we have string representations of the features corresponding to
  // those indices.
  // See sbn_maps.hpp for more description of these indexer representations.
  StringSetVector correct_representation_1(
      // The indexer representations for each of the possible virtual rootings.
      // For example, this first one is for rooting at the edge leading to leaf
      // 0, the second for rooting at leaf 1, etc.
      {{"00000|11111|01111", "10000|01111|00001", "00001|01110|00100",
        "00100|01010|00010"},
       {"00000|11111|01000", "01000|10111|00010", "00100|10001|00001",
        "00010|10101|00100"},
       {"00000|11111|00100", "10001|01010|00010", "01010|10001|00001",
        "00100|11011|01010"},
       {"00000|11111|00010", "00010|11101|01000", "00100|10001|00001",
        "01000|10101|00100"},
       {"00000|11111|00001", "00001|11110|01110", "10000|01110|00100",
        "00100|01010|00010"},
       {"00000|11111|01010", "10101|01010|00010", "00100|10001|00001",
        "01010|10101|00100"},
       {"00000|11111|01110", "00100|01010|00010", "10001|01110|00100",
        "01110|10001|00001"}});
  CHECK_EQ(
      inst.StringIndexerRepresentationOf(indexer_test_topology_1, out_of_sample_index),
      correct_representation_1);
  // See the "concepts" part of the online documentation to learn about PSP indexing.
  auto correct_psp_representation_1 =
      StringVectorVector({{"10000|01111", "10111|01000", "11011|00100", "11101|00010",
                           "11110|00001", "10101|01010", "10001|01110"},
                          {"", "", "", "", "", "01000|00010", "10000|00001"},
                          {"01110|00001", "10101|00010", "10001|01010", "10101|01000",
                           "10000|01110", "10001|00100", "01010|00100"}});
  CHECK_EQ(inst.psp_indexer_.StringRepresentationOf(indexer_test_topology_1),
           correct_psp_representation_1);
  // Same as above but for (((0,1),2),3,4);, or with internal nodes (((0,1)5,2)6,3,4)7;
  auto indexer_test_topology_2 = Node::OfParentIdVector({5, 5, 6, 7, 7, 6, 7});
  StringSetVector correct_representation_2(
      {{"00000|11111|01111", "10000|01111|00111", "00100|00011|00001",
        "01000|00111|00011"},
       {"00000|11111|01000", "01000|10111|00111", "00100|00011|00001",
        "10000|00111|00011"},
       {"00000|11111|00100", "00100|11011|00011", "11000|00011|00001",
        "00011|11000|01000"},
       {"00000|11111|00010", "00100|11000|01000", "00001|11100|00100",
        "00010|11101|00001"},
       {"00000|11111|00001", "00100|11000|01000", "00001|11110|00010",
        "00010|11100|00100"},
       {"00000|11111|00111", "00111|11000|01000", "00100|00011|00001",
        "11000|00111|00011"},
       {"00000|11111|00011", "00100|11000|01000", "11100|00011|00001",
        "00011|11100|00100"}});
  CHECK_EQ(
      inst.StringIndexerRepresentationOf(indexer_test_topology_2, out_of_sample_index),
      correct_representation_2);
  auto correct_psp_representation_2 =
      StringVectorVector({{"10000|01111", "10111|01000", "11011|00100", "11101|00010",
                           "11110|00001", "11000|00111", "11100|00011"},
                          {"", "", "", "", "", "10000|01000", "11000|00100"},
                          {"01000|00111", "10000|00111", "11000|00011", "11100|00001",
                           "11100|00010", "00100|00011", "00010|00001"}});
  CHECK_EQ(inst.psp_indexer_.StringRepresentationOf(indexer_test_topology_2),
           correct_psp_representation_2);

  // Test of RootedSBNMaps::IndexerRepresentationOf.
  // It's a little surprising to see this here in unrooted land, but these are actually
  // complementary tests to those found in rooted_sbn_instance.hpp, with a larger
  // subsplit support because we deroot the trees.
  // Topology is ((((0,1),2),3),4);, or with internal nodes ((((0,1)5,2)6,3)7,4)8;
  auto indexer_test_rooted_topology_1 =
      Node::OfParentIdVector({5, 5, 6, 7, 8, 6, 7, 8});
  auto correct_rooted_indexer_representation_1 =
      StringSet({"00000|11111|00001", "00001|11110|00010", "00010|11100|00100",
                 "00100|11000|01000"});
  CHECK_EQ(inst.StringIndexerRepresentationOf({RootedSBNMaps::IndexerRepresentationOf(
               inst.SBNSupport().Indexer(), indexer_test_rooted_topology_1,
               out_of_sample_index)})[0],
           correct_rooted_indexer_representation_1);
  // Topology is (((0,1),2),(3,4));, or with internal nodes (((0,1)5,2)6,(3,4)7)8;
  auto indexer_test_rooted_topology_2 =
      Node::OfParentIdVector({5, 5, 6, 7, 7, 6, 8, 8});
  auto correct_rooted_indexer_representation_2 =
      StringSet({"00000|11111|00011", "11100|00011|00001", "00011|11100|00100",
                 "00100|11000|01000"});
  CHECK_EQ(inst.StringIndexerRepresentationOf({RootedSBNMaps::IndexerRepresentationOf(
               inst.SBNSupport().Indexer(), indexer_test_rooted_topology_2,
               out_of_sample_index)})[0],
           correct_rooted_indexer_representation_2);
}

TEST_CASE("UnrootedSBNInstance: likelihood and gradient") {
  UnrootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/hello.nwk");
  inst.ReadFastaFile("data/hello.fasta");
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  inst.PrepareForPhyloLikelihood(simple_specification, 2);
  for (auto ll : inst.LogLikelihoods()) {
    CHECK_LT(fabs(ll - -84.852358), 0.000001);
  }
  inst.ReadNexusFile("data/DS1.subsampled_10.t");
  inst.ReadFastaFile("data/DS1.fasta");
  std::vector<BeagleFlags> vector_flag_options{BEAGLE_FLAG_VECTOR_NONE,
                                               BEAGLE_FLAG_VECTOR_SSE};
  std::vector<bool> tip_state_options{false, true};
  for (const auto vector_flag : vector_flag_options) {
    for (const auto tip_state_option : tip_state_options) {
      inst.PrepareForPhyloLikelihood(simple_specification, 2, {vector_flag},
                                     tip_state_option);
      auto likelihoods = inst.LogLikelihoods();
      std::vector<double> pybeagle_likelihoods(
          {-14582.995273982739, -6911.294207416366, -6916.880235529542,
           -6904.016888831189, -6915.055570693576, -6915.50496696512,
           -6910.958836661867, -6909.02639968063, -6912.967861935749,
           -6910.7871105783515});
      for (size_t i = 0; i < likelihoods.size(); i++) {
        CHECK_LT(fabs(likelihoods[i] - pybeagle_likelihoods[i]), 0.00011);
      }

      auto gradients = inst.PhyloGradients();
      // Test the log likelihoods.
      for (size_t i = 0; i < likelihoods.size(); i++) {
        CHECK_LT(fabs(gradients[i].log_likelihood_ - pybeagle_likelihoods[i]), 0.00011);
      }
      // Test the gradients for the last tree.
      auto last = gradients.back();
      std::sort(last.gradient_["branch_lengths"].begin(),
                last.gradient_["branch_lengths"].end());
      // Zeros are for the root and one of the descendants of the root.
      std::vector<double> physher_gradients = {
          -904.18956, -607.70500, -562.36274, -553.63315, -542.26058, -539.64210,
          -463.36511, -445.32555, -414.27197, -412.84218, -399.15359, -342.68038,
          -306.23644, -277.05392, -258.73681, -175.07391, -171.59627, -168.57646,
          -150.57623, -145.38176, -115.15798, -94.86412,  -83.02880,  -80.09165,
          -69.00574,  -51.93337,  0.00000,    0.00000,    16.17497,   20.47784,
          58.06984,   131.18998,  137.10799,  225.73617,  233.92172,  253.49785,
          255.52967,  259.90378,  394.00504,  394.96619,  396.98933,  429.83873,
          450.71566,  462.75827,  471.57364,  472.83161,  514.59289,  650.72575,
          888.87834,  913.96566,  927.14730,  959.10746,  2296.55028};
      for (size_t i = 0; i < last.gradient_["branch_lengths"].size(); i++) {
        CHECK_LT(fabs(last.gradient_["branch_lengths"][i] - physher_gradients[i]),
                 0.0001);
      }

      // Test rescaling
      inst.SetRescaling(true);
      auto likelihoods_rescaling = inst.LogLikelihoods();
      // Likelihoods from LogLikelihoods()
      for (size_t i = 0; i < likelihoods_rescaling.size(); i++) {
        CHECK_LT(fabs(likelihoods_rescaling[i] - pybeagle_likelihoods[i]), 0.00011);
      }
      // Likelihoods from BranchGradients()
      inst.PrepareForPhyloLikelihood(simple_specification, 1, {}, tip_state_option);
      auto gradients_rescaling = inst.PhyloGradients();
      for (size_t i = 0; i < gradients_rescaling.size(); i++) {
        CHECK_LT(fabs(gradients_rescaling[i].log_likelihood_ - pybeagle_likelihoods[i]),
                 0.00011);
      }
      // Gradients
      auto last_rescaling = gradients_rescaling.back();
      auto branch_lengths_gradient = last_rescaling.gradient_["branch_lengths"];
      std::sort(branch_lengths_gradient.begin(), branch_lengths_gradient.end());
      for (size_t i = 0; i < branch_lengths_gradient.size(); i++) {
        CHECK_LT(fabs(branch_lengths_gradient[i] - physher_gradients[i]), 0.0001);
      }
    }
  }
}

TEST_CASE("UnrootedSBNInstance: likelihood and gradient with Weibull") {
  UnrootedSBNInstance inst("charlie");
  PhyloModelSpecification simple_specification{"JC69", "weibull+4", "strict"};
  inst.ReadNexusFile("data/DS1.subsampled_10.t");
  inst.ReadFastaFile("data/DS1.fasta");

  std::vector<double> physher_likelihoods(
      {-9456.1201098061, -6624.4110704332, -6623.4474776131, -6617.25658038029,
       -6627.5385571548, -6621.6155048722, -6622.3314942713, -6618.7695717585,
       -6616.3837517370, -6623.8295828648});
  // First element of each gradient
  std::vector<double> physher_gradients_bl0(
      {-126.890527, 157.251275, 138.202510, -180.311856, 417.562897, -796.450894,
       -173.744375, -70.693513, 699.190754, -723.034349});
  std::vector<BeagleFlags> vector_flag_options{BEAGLE_FLAG_VECTOR_NONE,
                                               BEAGLE_FLAG_VECTOR_SSE};
  std::vector<bool> tip_state_options{false, true};
  for (const auto vector_flag : vector_flag_options) {
    for (const auto tip_state_option : tip_state_options) {
      inst.PrepareForPhyloLikelihood(simple_specification, 2, {vector_flag},
                                     tip_state_option);
      auto param_block_map = inst.GetPhyloModelParamBlockMap();
      param_block_map.at(WeibullSiteModel::shape_key_).setConstant(0.1);
      auto likelihoods = inst.LogLikelihoods();
      for (size_t i = 0; i < likelihoods.size(); i++) {
        CHECK_LT(fabs(likelihoods[i] - physher_likelihoods[i]), 0.00011);
      }

      auto gradients = inst.PhyloGradients();
      for (size_t i = 0; i < gradients.size(); i++) {
        CHECK_LT(fabs(gradients[i].gradient_["branch_lengths"][0] -
                      physher_gradients_bl0[i]),
                 0.00011);
      }

      // Test rescaling
      inst.SetRescaling(true);
      auto likelihoods_rescaling = inst.LogLikelihoods();
      // Likelihoods from LogLikelihoods()
      for (size_t i = 0; i < likelihoods_rescaling.size(); i++) {
        CHECK_LT(fabs(likelihoods_rescaling[i] - physher_likelihoods[i]), 0.00011);
      }

      auto gradients_rescaling = inst.PhyloGradients();
      for (size_t i = 0; i < gradients.size(); i++) {
        CHECK_LT(fabs(gradients_rescaling[i].gradient_["branch_lengths"][0] -
                      physher_gradients_bl0[i]),
                 0.00011);
      }
    }
  }
}

TEST_CASE("UnrootedSBNInstance: SBN training") {
  UnrootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/DS1.100_topologies.nwk");
  inst.ProcessLoadedTrees();
  // These "Expected" functions are defined in sbn_probability.hpp.
  const auto expected_SA = ExpectedSAVector();
  inst.TrainSimpleAverage();
  CheckVectorXdEquality(inst.CalculateSBNProbabilities(), expected_SA, 1e-12);
  // Expected EM vectors with alpha = 0.
  const auto [expected_EM_0_1, expected_EM_0_23] = ExpectedEMVectorsAlpha0();
  // 1 iteration of EM with alpha = 0.
  inst.TrainExpectationMaximization(0., 1);
  CheckVectorXdEquality(inst.CalculateSBNProbabilities(), expected_EM_0_1, 1e-12);
  // 23 iterations of EM with alpha = 0.
  inst.TrainExpectationMaximization(0., 23);
  CheckVectorXdEquality(inst.CalculateSBNProbabilities(), expected_EM_0_23, 1e-12);
  // 100 iteration of EM with alpha = 0.5.
  const auto expected_EM_05_100 = ExpectedEMVectorAlpha05();
  inst.TrainExpectationMaximization(0.5, 100);
  CheckVectorXdEquality(inst.CalculateSBNProbabilities(), expected_EM_05_100, 1e-5);
}

TEST_CASE("UnrootedSBNInstance: tree sampling") {
  UnrootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon_unrooted.nwk");
  inst.ProcessLoadedTrees();
  inst.TrainSimpleAverage();
  // Count the frequencies of rooted trees in a file.
  size_t rooted_tree_count_from_file = 0;
  RootedIndexerRepresentationSizeDict counter_from_file(0);
  for (const auto &indexer_representation : inst.MakeIndexerRepresentations()) {
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(counter_from_file,
                                                                indexer_representation);
    rooted_tree_count_from_file += indexer_representation.size();
  }
  // Count the frequencies of trees when we sample after training with
  // SimpleAverage.
  size_t sampled_tree_count = 1'000'000;
  RootedIndexerRepresentationSizeDict counter_from_sampling(0);
  ProgressBar progress_bar(sampled_tree_count / 1000);
  for (size_t sample_idx = 0; sample_idx < sampled_tree_count; ++sample_idx) {
    const auto rooted_topology = inst.SampleTopology(true);
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(
        counter_from_sampling,
        RootedSBNMaps::IndexerRepresentationOf(inst.SBNSupport().Indexer(),
                                               rooted_topology, out_of_sample_index));
    if (sample_idx % 1000 == 0) {
      ++progress_bar;
      progress_bar.display();
    }
  }
  // These should be equal in the limit when we're training with SA.
  for (const auto &[key, _] : counter_from_file) {
    std::ignore = _;
    double observed =
        static_cast<double>(counter_from_sampling.at(key)) / sampled_tree_count;
    double expected =
        static_cast<double>(counter_from_file.at(key)) / rooted_tree_count_from_file;
    CHECK_LT(fabs(observed - expected), 5e-3);
  }
  progress_bar.done();
}

TEST_CASE("UnrootedSBNInstance: gradient of log q_{phi}(tau) WRT phi") {
  UnrootedSBNInstance inst("charlie");
  // File gradient_test.t contains two trees:
  // ((0,1), 2, (3,4)) and
  // ((0,1), (2,3), 4).
  inst.ReadNexusFile("data/gradient_test.t");
  inst.ProcessLoadedTrees();

  // The number of rootsplits across all of the input trees.
  size_t num_rootsplits = 8;
  // Manual enumeration shows that there are 31 PCSP's.
  size_t num_pcsp = inst.sbn_parameters_.size() - num_rootsplits;

  // Test for K = 1 tree.
  size_t K = 1;
  inst.tree_collection_.trees_.clear();
  // Generate a tree,
  // \tau = ((0,1),(2,3),4) with internal node labels ((0,1)5,(2,3)6,4)7.
  std::vector<size_t> tau_indices = {5, 5, 6, 6, 7, 7, 7};
  auto tau = UnrootedTree::OfParentIdVector(tau_indices);
  inst.tree_collection_.trees_.push_back(tau);

  // Initialize sbn_parameters_ to 0's and normalize, which is going to give a uniform
  // distribution for rootsplits and PCSP distributions.
  inst.sbn_parameters_.setZero();
  EigenVectorXd normalized_sbn_parameters_in_log = inst.sbn_parameters_;
  inst.ProbabilityNormalizeSBNParametersInLog(normalized_sbn_parameters_in_log);
  // Because this is a uniform distribution, each rootsplit \rho has P(\rho) = 1/8.
  //
  // We're going to start by computing the rootsplit gradient.
  // There are 7 possible rootings of \tau.
  // For example consider rooting on the 014|23 split, yielding the following subsplits:
  // 014|23, 2|3, 01|4, 0|1.
  // Each of the child subsplits are the only possible subsplit,
  // except for the root where it has probability 1/8. Hence, the probability
  // for this tree is 1/8 x 1 x 1 x 1 = 1/8.
  // Now, consider rooting on the 0|1234 split, yielding the following subsplits:
  // 0|1234, 1|234, 23|4, 2|3.
  // The probability for this tree is 1/8 x 1 x 1/2 x 1 = 1/16, where the 1/2 comes from
  // the fact that we can have 23|4 or 2|34.
  //
  // Each of the remaining 5 trees has the same probability: the product of
  // 1/8 for the rootsplit and 1/2 for one of the subsplit resolutions of 234.
  // One can see this because the only way for there not to be ambiguity in the
  // resolution of the splitting of 234 is for one to take 014|23 as the rootsplit.
  //
  // Hence, q(\tau) = 6 x 1/16 + 1 x 1/8 = 8/16 = 0.5.
  // Note that there are a total of 8 rootsplits; 7 are possible rootsplits of
  // the sampled tree \tau but one rootsplit,
  // 014|23 is not observed rooting of \tau and hence,
  // the gradient for 014|23 is simply -P(014|23) = -1/8.
  //
  // The gradient with respect to each of the 7 rootsplits is given by
  // P(\tau_{\rho})/q(\tau) - P(\rho) via eq:rootsplitGrad,
  // which is equal to
  // (1/8) / (0.5) - 1/8 = 1/8 for the tree with \rho = 34|125 and
  // (1/16) / (0.5) - 1/8 = 0 for 6 remaining trees.

  EigenVectorXd expected_grad_rootsplit(8);
  expected_grad_rootsplit << -1. / 8, 0, 0, 0, 0, 0, 0, 1. / 8;
  auto indexer_representations = inst.MakeIndexerRepresentations();
  EigenVectorXd grad_log_q = inst.GradientOfLogQ(normalized_sbn_parameters_in_log,
                                                 indexer_representations.at(0));
  EigenVectorXd realized_grad_rootsplit = grad_log_q.segment(0, 8);
  // Sort them and compare against sorted version of
  // realized_grad_rootsplit[0:7].
  std::sort(realized_grad_rootsplit.begin(), realized_grad_rootsplit.end());
  CheckVectorXdEquality(realized_grad_rootsplit, expected_grad_rootsplit, 1e-8);

  // Manual enumeration shows that the entries corresponding to PCSP should have
  // 6 entries with -1/16 and 6 entries with 1/16 and the rest with 0's.
  // For example, consider the tree ((0,1),(2,3),4), which has the following subsplits:
  // 0123|4, 01|23, 0|1, 2|3.
  // Note the subsplit s = 01|23 is one of two choices for
  // the parent subsplit t = 0123|4,
  // since 0123|4 can also be split into s' = 012|3.
  // Let \rho = 0123|4, the gradient for 01|23 is given by:
  // (1/q(\tau)) P(\tau_{\rho}) * (1 - P(01|23 | 0123|4))
  // = 2 * (1/16) * (1-0.5) = 1/16.
  // The gradient for s' = 012|3 is,
  // (1/q(\tau)) P(\tau_{\rho}) * -P(012|3 | 0123|4)
  // = 2 * (1/16) * -0.5 = -1/16.

  // The gradient for the following PCSP are 1/16 as above.
  // 014|3 | 0134|2
  // 014|2 | 0124|3
  // 01|23 | 0123|4
  // 23|4  | 01|234
  // 23|4  | 1|234
  // 23|4  | 0|234
  // And each of these have an alternate subsplit s' that gets a gradient of -1/16.
  // Each of the other PCSP gradients are 0 either because its parent support
  // never appears in the tree or it represents the only child subsplit.
  EigenVectorXd expected_grad_pcsp = EigenVectorXd::Zero(num_pcsp);
  expected_grad_pcsp.segment(0, 6).setConstant(-1. / 16);
  expected_grad_pcsp.segment(num_pcsp - 6, 6).setConstant(1. / 16);
  EigenVectorXd realized_grad_pcsp = grad_log_q.tail(num_pcsp);
  std::sort(realized_grad_pcsp.begin(), realized_grad_pcsp.end());
  CheckVectorXdEquality(realized_grad_pcsp, expected_grad_pcsp, 1e-8);

  // We'll now change the SBN parameters and check the gradient there.
  // If we root at 0123|4, then the only choice we have is between the following s and
  // s' as described above.
  // The PCSP s|t =  (01|23) | (0123|4) corresponds to 00001|11110|00110.
  // The PCSP s'|t = (012|3) | (0123|4) corresponds to 00001|11110|00010.
  Bitset s("000011111000110");
  Bitset s_prime("000011111000010");
  size_t s_idx = inst.SBNSupport().IndexerAt(s);
  size_t s_prime_idx = inst.SBNSupport().IndexerAt(s_prime);
  inst.sbn_parameters_.setZero();
  inst.sbn_parameters_(s_idx) = 1;
  inst.sbn_parameters_(s_prime_idx) = -1;
  normalized_sbn_parameters_in_log = inst.sbn_parameters_;
  inst.ProbabilityNormalizeSBNParametersInLog(normalized_sbn_parameters_in_log);

  // These changes to normalized_sbn_parameters_in_log will change q(\tau) as well as
  // P(\tau_{\rho}) for \rho = 0123|4. First,
  // P(\tau_{\rho}) = 1/8 * exp(1)/(exp(1) + exp(-1)) = 0.1100996.
  double p_tau_rho = (1. / 8) * exp(normalized_sbn_parameters_in_log[s_idx]);
  // For q(\tau), we will just compute using the already tested function:
  double q_tau = inst.CalculateSBNProbabilities()(0);
  // The gradient for s|t is given by,
  // (1/q(\tau)) x P(\tau_{\rho}) x (1 - P(s|t))
  double expected_grad_at_s =
      (1. / q_tau) * p_tau_rho * (1 - exp(normalized_sbn_parameters_in_log[s_idx]));
  // And the gradient for s'|t is given by,
  // (1/q(\tau)) x P(\tau_{\rho}) x (-P(s|t))
  double expected_grad_at_s_prime =
      (1. / q_tau) * p_tau_rho * -exp(normalized_sbn_parameters_in_log[s_prime_idx]);
  // We're setting normalized_sbn_parameters_in_log to NaN as we would in a normal
  // application of GradientOfLogQ.
  normalized_sbn_parameters_in_log.setConstant(DOUBLE_NAN);
  grad_log_q = inst.GradientOfLogQ(normalized_sbn_parameters_in_log,
                                   indexer_representations.at(0));
  CHECK_LT(fabs(expected_grad_at_s - grad_log_q(s_idx)), 1e-8);
  CHECK_LT(fabs(expected_grad_at_s_prime - grad_log_q(s_prime_idx)), 1e-8);

  // Now we test the gradient by doing the calculation by hand.
  K = 4;
  inst.SampleTrees(K);
  // Make up some numbers for log_f.
  EigenVectorXd log_f(K);
  log_f << -83, -75, -80, -79;
  // log_F = -74.97493
  double log_F = NumericalUtils::LogSum(log_f);
  double elbo = log_F - log(K);
  // 0.0003271564 0.9752395946 0.0065711127 0.0178621362
  EigenVectorXd tilde_w = (log_f.array() - log_F).exp();
  // -76.36155 -77.33646 -76.36779 -76.37908
  EigenVectorXd multiplicative_factors = (elbo - tilde_w.array());

  EigenVectorXd expected_nabla(inst.sbn_parameters_.size());
  expected_nabla.setZero();
  // We now have some confidence in GradientOfLogQ(), so we just use it.
  auto indexer_reps = inst.MakeIndexerRepresentations();
  normalized_sbn_parameters_in_log.setConstant(DOUBLE_NAN);
  for (size_t k = 0; k < K; k++) {
    grad_log_q =
        multiplicative_factors(k) *
        inst.GradientOfLogQ(normalized_sbn_parameters_in_log, indexer_reps.at(k))
            .array();
    expected_nabla += grad_log_q;
  }
  bool use_vimco = false;
  EigenVectorXd realized_nabla = inst.TopologyGradients(log_f, use_vimco);
  CheckVectorXdEquality(realized_nabla, expected_nabla, 1e-8);

  // Test for VIMCO gradient estimator.
  EigenVectorXd vimco_multiplicative_factors(K);
  vimco_multiplicative_factors << -0.04742748, 2.59553236, -0.01779887, -0.01278592;
  expected_nabla.setZero();
  normalized_sbn_parameters_in_log.setConstant(DOUBLE_NAN);
  for (size_t k = 0; k < K; k++) {
    grad_log_q =
        vimco_multiplicative_factors(k) *
        inst.GradientOfLogQ(normalized_sbn_parameters_in_log, indexer_reps.at(k))
            .array();
    expected_nabla += grad_log_q;
  }
  use_vimco = true;
  realized_nabla = inst.TopologyGradients(log_f, use_vimco);
  CheckVectorXdEquality(realized_nabla, expected_nabla, 1e-8);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
