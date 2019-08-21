// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIBSBN_HPP_
#define SRC_LIBSBN_HPP_

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include "alignment.hpp"
#include "beagle.hpp"
#include "driver.hpp"
#include "prettyprint.hpp"
#include "psp_indexer.hpp"
#include "sbn_maps.hpp"
#include "sugar.hpp"
#include "tree.hpp"

typedef std::unordered_map<std::string, float> StringFloatMap;
typedef std::unordered_map<std::string, uint32_t> StringUInt32Map;
typedef std::unordered_map<std::string, std::pair<uint32_t, uint32_t>>
    StringUInt32PairMap;
typedef std::unordered_map<uint32_t, std::string> UInt32StringMap;
typedef std::unordered_map<std::string,
                           std::unordered_map<std::string, uint32_t>>
    StringPCSSMap;
typedef std::vector<uint32_t> PCSSIndexVector;

// Turn a <Key, T> map into a <std::string, T> map for any Key type that has a
// ToString method.
template <class Key, class T>
std::unordered_map<std::string, T> StringifyMap(std::unordered_map<Key, T> m) {
  std::unordered_map<std::string, T> m_str;
  for (const auto &iter : m) {
    m_str[iter.first.ToString()] = iter.second;
  }
  return m_str;
}

StringPCSSMap StringPCSSMapOf(PCSSDict d);

struct SBNInstance {
  std::string name_;
  // Trees get loaded in from a file or sampled from SBNs.
  TreeCollection tree_collection_;
  Alignment alignment_;
  // Beagly bits.
  CharIntMap symbol_table_;
  std::vector<beagle::BeagleInstance> beagle_instances_;
  size_t beagle_leaf_count_;
  size_t beagle_site_count_;
  // A vector that contains all of the SBN-related probabilities.
  std::vector<double> sbn_parameters_;
  // A map that indexes these probabilities: rootsplits are at the beginning,
  // and PCSS bitsets are at the end.
  // The collection of rootsplits, with the same indexing as in the indexer_.
  BitsetVector rootsplits_;
  // The first index after the rootsplit block in sbn_parameters_.
  size_t rootsplit_index_end_;
  BitsetUInt32Map indexer_;
  // A map going from the index of a PCSS to its child.
  UInt32BitsetMap index_to_child_;
  // A map going from a parent subsplit to the range of indices in
  // sbn_parameters_ with its children.
  BitsetUInt32PairMap parent_to_range_;
  // Random bits.
  static std::random_device random_device_;
  static std::mt19937 random_generator_;
  bool rescaling_;

  // ** Initialization, destruction, and status
  explicit SBNInstance(const std::string &name)
      : name_(name),
        symbol_table_(beagle::GetSymbolTable()),
        beagle_leaf_count_(0),
        beagle_site_count_(0),
        rescaling_{false} {}

  ~SBNInstance() { FinalizeBeagleInstances(); }

  // Release memory from BEAGLE instances and discard them.
  void FinalizeBeagleInstances();

  void SetRescaling(bool use_rescaling) { rescaling_ = use_rescaling; }

  size_t TreeCount() const { return tree_collection_.TreeCount(); }
  void PrintStatus();

  // ** Building SBN-related items
  // Define "SBN maps" to be the collection of maps associated with the
  // SBNInstance, such as indexer_, index_to_child_, parent_to_range_, and
  // rootsplits_.

  // Use the loaded trees to get the SBN maps and prepare the sbn_parameters_
  // vector.
  void ProcessLoadedTrees();
  void CheckSBNMapsAvailable();
  void PrintSupports();

  // Sample an integer index in [range.first, range.second) according to
  // sbn_parameters_.
  uint32_t SampleIndex(std::pair<uint32_t, uint32_t> range) const;

  // This function samples a tree by first sampling the rootsplit, and then
  // calling the recursive form of SampleTopology (see next function)..
  Node::NodePtr SampleTopology() const;

  // The input to this function is a parent subsplit (of length 2n).
  Node::NodePtr SampleTopology(const Bitset &parent_subsplit) const;

  void SampleTrees(size_t count);

  // Get indexer representations of the trees in tree_collection_.
  // See the header documentation of IndexerRepresentationOf for an explanation
  // of what these are.
  std::vector<IndexerRepresentation> GetIndexerRepresentations() const;

  // Get the indexer, but reversed and with bitsets appropriately converted to
  // strings.
  StringVector StringReversedIndexer() const;

  // Turn an IndexerRepresentation into a string representation of the underying
  // bitsets. This is really just so that we can make a test of indexer
  // representations.
  std::pair<StringSet, StringSetVector> StringIndexerRepresentationOf(
      IndexerRepresentation indexer_representation) const;

  // ** I/O

  // Return indexer_ and parent_to_range_ converted into string-keyed maps.
  std::tuple<StringUInt32Map, StringUInt32PairMap> GetIndexers() const;

  // This function is really just for testing-- it recomputes counters scratch.
  std::pair<StringUInt32Map, StringPCSSMap> SplitCounters() const;

  void ReadNewickFile(std::string fname);
  void ReadNexusFile(std::string fname);
  void ReadFastaFile(std::string fname);

  // ** Phylogenetic likelihood

  void CheckSequencesAndTreesLoaded() const;

  // Make sure that BEAGLE's idea of the sequence count and length match what we
  // have loaded.
  void CheckBeagleDimensions() const;

  // Make the specified number of BEAGLE instances are appropriate for the
  // loaded data. Likelihood calculation will be run in parallel across the
  // specified number of instances.
  void MakeBeagleInstances(int instance_count);

  std::vector<double> LogLikelihoods() const;
  // For each loaded tree, returns a pair of (likelihood, gradient).
  std::vector<std::pair<double, std::vector<double>>> BranchGradients() const;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("libsbn") {
  SBNInstance inst("charlie");
  inst.ReadNewickFile("data/hello.nwk");
  inst.ReadFastaFile("data/hello.fasta");
  inst.MakeBeagleInstances(2);
  for (auto ll : inst.LogLikelihoods()) {
    CHECK_LT(fabs(ll - -84.852358), 0.000001);
  }
  // Reading one file after another checks that we've cleared out state.
  inst.ReadNewickFile("data/five_taxon.nwk");
  inst.ProcessLoadedTrees();
  // (2,(1,3),(0,4));, or with internal nodes (2,(1,3)5,(0,4)6)7
  auto indexer_test_topology_1 = Node::OfParentIdVector({6, 5, 7, 5, 6, 7, 7});
  std::pair<StringSet, StringSetVector> correct_representation_1(
      // The rootsplits.
      {"01110", "01000", "01010", "01111", "00010", "00100", "00001"},
      // The PCSSs for each of the possible virtual rootings.
      // For example, this first one is for rooting at the edge leading to leaf
      // 0.
      {{"10000|01111|00001", "00001|01110|00100", "00100|01010|00010"},
       {"01000|10111|00010", "00100|10001|00001", "00010|10101|00100"},
       {"10001|01010|00010", "01010|10001|00001", "00100|11011|01010"},
       {"00010|11101|01000", "00100|10001|00001", "01000|10101|00100"},
       {"00001|11110|01110", "10000|01110|00100", "00100|01010|00010"},
       {"10101|01010|00010", "00100|10001|00001", "01010|10101|00100"},
       {"00100|01010|00010", "10001|01110|00100", "01110|10001|00001"}});
  CHECK_EQ(inst.StringIndexerRepresentationOf(
               IndexerRepresentationOf(inst.indexer_, indexer_test_topology_1)),
           correct_representation_1);
  // (((0,1),2),3,4);, or with internal nodes (((0,1)5,2)6,3,4)7;
  auto indexer_test_topology_2 = Node::OfParentIdVector({5, 5, 6, 7, 7, 6, 7});
  std::pair<StringSet, StringSetVector> correct_representation_2(
      {"01000", "01111", "00011", "00010", "00111", "00100", "00001"},
      {{"10000|01111|00111", "00100|00011|00001", "01000|00111|00011"},
       {"01000|10111|00111", "00100|00011|00001", "10000|00111|00011"},
       {"00100|11011|00011", "11000|00011|00001", "00011|11000|01000"},
       {"00100|11000|01000", "00001|11100|00100", "00010|11101|00001"},
       {"00100|11000|01000", "00001|11110|00010", "00010|11100|00100"},
       {"00111|11000|01000", "00100|00011|00001", "11000|00111|00011"},
       {"00100|11000|01000", "11100|00011|00001", "00011|11100|00100"}});
  CHECK_EQ(inst.StringIndexerRepresentationOf(
               IndexerRepresentationOf(inst.indexer_, indexer_test_topology_2)),
           correct_representation_2);
  inst.SampleTrees(2);
  inst.GetIndexerRepresentations();

  std::cout << PSPRepresentationOf(inst.indexer_, indexer_test_topology_1)
            << std::endl;

  inst.ReadNexusFile("data/DS1.subsampled_10.t");
  inst.ReadFastaFile("data/DS1.fasta");
  inst.MakeBeagleInstances(2);
  auto likelihoods = inst.LogLikelihoods();
  std::vector<double> pybeagle_likelihoods(
      {-14582.995273982739, -6911.294207416366, -6916.880235529542,
       -6904.016888831189, -6915.055570693576, -6915.50496696512,
       -6910.958836661867, -6909.02639968063, -6912.967861935749,
       -6910.7871105783515});
  for (size_t i = 0; i < likelihoods.size(); i++) {
    CHECK_LT(fabs(likelihoods[i] - pybeagle_likelihoods[i]), 0.00011);
  }

  auto gradients = inst.BranchGradients();
  // Test the log likelihoods.
  for (size_t i = 0; i < likelihoods.size(); i++) {
    CHECK_LT(fabs(gradients[i].first - pybeagle_likelihoods[i]), 0.00011);
  }
  // Test the gradients for the last tree.
  auto last = gradients.back();
  std::sort(last.second.begin(), last.second.end());
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
  for (size_t i = 0; i < last.second.size(); i++) {
    CHECK_LT(fabs(last.second[i] - physher_gradients[i]), 0.0001);
  }

  // Test rescaling
  inst.SetRescaling(true);
  auto likelihoods_rescaling = inst.LogLikelihoods();
  // Likelihoods from LogLikelihoods()
  for (size_t i = 0; i < likelihoods_rescaling.size(); i++) {
    CHECK_LT(fabs(likelihoods_rescaling[i] - pybeagle_likelihoods[i]), 0.00011);
  }
  // Likelihoods from BranchGradients()
  inst.MakeBeagleInstances(1);
  auto gradients_rescaling = inst.BranchGradients();
  for (size_t i = 0; i < gradients_rescaling.size(); i++) {
    CHECK_LT(fabs(gradients_rescaling[i].first - pybeagle_likelihoods[i]),
             0.00011);
  }
  // Gradients
  auto last_rescaling = gradients_rescaling.back();
  std::sort(last_rescaling.second.begin(), last_rescaling.second.end());
  for (size_t i = 0; i < last_rescaling.second.size(); i++) {
    CHECK_LT(fabs(last_rescaling.second[i] - physher_gradients[i]), 0.0001);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_LIBSBN_HPP_
