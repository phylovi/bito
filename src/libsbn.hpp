// Copyright 2019 Matsen group.
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
#include "build.hpp"
#include "driver.hpp"
#include "prettyprint.hpp"
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

template <typename T>
StringUInt32Map StringUInt32MapOf(T m) {
  StringUInt32Map m_str;
  for (const auto &iter : m) {
    m_str[iter.first.ToString()] = iter.second;
  }
  return m_str;
}

StringPCSSMap StringPCSSMapOf(PCSSDict d) {
  StringPCSSMap d_str;
  for (const auto &iter : d) {
    d_str[iter.first.ToString()] = StringUInt32MapOf(iter.second);
  }
  return d_str;
}

struct SBNInstance {
  std::string name_;
  // Things that get loaded in or sampled from SBNs.
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

  // Finalize means to release memory.
  void FinalizeBeagleInstances() {
    for (const auto &beagle_instance : beagle_instances_) {
      assert(beagleFinalizeInstance(beagle_instance) == 0);
    }
    beagle_instances_.clear();
    beagle_leaf_count_ = 0;
    beagle_site_count_ = 0;
  }

  void SetRescaling(bool use_rescaling) { rescaling_ = use_rescaling; }

  size_t TreeCount() const { return tree_collection_.TreeCount(); }
  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (tree_collection_.TreeCount()) {
      std::cout << TreeCount() << " unique tree topologies loaded on "
                << tree_collection_.TaxonCount() << " leaves.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
    std::cout << alignment_.Data().size() << " sequences loaded.\n";
  }

  // ** Building SBN-related items

  void ProcessLoadedTrees() {
    uint32_t index = 0;
    auto counter = tree_collection_.TopologyCounter();
    // See above for the definitions of these members.
    sbn_parameters_.clear();
    rootsplits_.clear();
    indexer_.clear();
    index_to_child_.clear();
    parent_to_range_.clear();
    // Start by adding the rootsplits.
    for (const auto &iter : RootsplitCounterOf(counter)) {
      assert(indexer_.insert({iter.first, index}).second);
      rootsplits_.push_back(iter.first);
      index++;
    }
    rootsplit_index_end_ = index;
    // Now add the PCSSs.
    for (const auto &iter : PCSSCounterOf(counter)) {
      const auto &parent = iter.first;
      const auto &child_counter = iter.second;
      assert(parent_to_range_
                 .insert({parent, {index, index + child_counter.size()}})
                 .second);
      for (const auto &child_iter : child_counter) {
        const auto &child = child_iter.first;
        assert(indexer_.insert({parent + child, index}).second);
        assert(index_to_child_
                   .insert({index, Bitset::ChildSubsplit(parent, child)})
                   .second);
        index++;
      }
    }
    sbn_parameters_ = std::vector<double>(index, 1.);
  }

  void PrintSupports() {
    std::vector<std::string> to_print(indexer_.size());
    for (const auto &iter : indexer_) {
      if (iter.second < rootsplit_index_end_) {
        to_print[iter.second] = iter.first.ToString();
      } else {
        to_print[iter.second] = iter.first.PCSSToString();
      }
    }
    for (size_t i = 0; i < to_print.size(); i++) {
      std::cout << i << "\t" << to_print[i] << std::endl;
    }
  }

  // Sample an integer index in [range.first, range.second) according to
  // sbn_parameters_.
  uint32_t SampleIndex(std::pair<uint32_t, uint32_t> range) const {
    assert(range.first < range.second);
    assert(range.second <= sbn_parameters_.size());
    std::discrete_distribution<> distribution(
        sbn_parameters_.begin() + range.first,
        sbn_parameters_.begin() + range.second);
    // We have to add on range.first because we have taken a slice of the full
    // array, and the sampler treats the beginning of this slice as zero.
    auto result =
        range.first + static_cast<uint32_t>(distribution(random_generator_));
    assert(result < range.second);
    return result;
  }

  // This function samples a tree by first sampling the rootsplit, and then
  // calling the recursive form of SampleTopology.
  Node::NodePtr SampleTopology() const {
    // Start by sampling a rootsplit.
    uint32_t rootsplit_index =
        SampleIndex(std::pair<uint32_t, uint32_t>(0, rootsplit_index_end_));
    const Bitset &rootsplit = rootsplits_.at(rootsplit_index);
    // The addition below turns the rootsplit into a subsplit.
    auto topology = SampleTopology(rootsplit + ~rootsplit)->Deroot();
    topology->Reindex();
    return topology;
  }

  // The input to this function is a parent subsplit (of length 2n).
  Node::NodePtr SampleTopology(const Bitset &parent_subsplit) const {
    auto process_subsplit = [this](const Bitset &parent) {
      auto singleton_option = parent.SplitChunk(1).SingletonOption();
      if (singleton_option) {
        return Node::Leaf(*singleton_option);
      }  // else
      auto child_index = SampleIndex(parent_to_range_.at(parent));
      return SampleTopology(index_to_child_.at(child_index));
    };
    return Node::Join(process_subsplit(parent_subsplit),
                      process_subsplit(parent_subsplit.RotateSubsplit()));
  }

  void CheckSBNMapsAvailable() {
    if (!indexer_.size() || !index_to_child_.size() ||
        !parent_to_range_.size() || !rootsplits_.size()) {
      std::cout << "Please call ProcessLoadedTrees to prepare your SBN maps.\n";
      abort();
    }
  }

  void SampleTrees(size_t count) {
    CheckSBNMapsAvailable();
    auto leaf_count = rootsplits_[0].size();
    // 2n-2 because trees are unrooted.
    auto edge_count = 2 * static_cast<int>(leaf_count) - 2;
    tree_collection_.trees_.clear();
    for (size_t i = 0; i < count; i++) {
      std::vector<double> branch_lengths(static_cast<size_t>(edge_count));
      tree_collection_.trees_.emplace_back(
          Tree(SampleTopology(), std::move(branch_lengths)));
    }
  }

  std::vector<IndexerRepresentation> GetIndexerRepresentations() {
    std::vector<IndexerRepresentation> representations;
    representations.reserve(tree_collection_.trees_.size());
    for (const auto &tree : tree_collection_.trees_) {
      IndexerRepresentationOf(indexer_, tree.Topology());
      representations.push_back(
          IndexerRepresentationOf(indexer_, tree.Topology()));
    }
    return representations;
  }

  // Get the indexer, but reversed and with bitsets appropriately converted to
  // strings.
  StringVector StringReversedIndexer() {
    std::vector<std::string> reversed_indexer(indexer_.size());
    for (const auto &iter : indexer_) {
      if (iter.second < rootsplit_index_end_) {
        reversed_indexer[iter.second] = iter.first.ToString();
      } else {
        reversed_indexer[iter.second] = iter.first.PCSSToString();
      }
    }
    return reversed_indexer;
  }

  // Turn an IndexerRepresentation into a string representation of the underying
  // bitsets.
  std::pair<StringSet, StringSetVector> StringIndexerRepresentationOf(
      IndexerRepresentation indexer_representation) {
    auto reversed_indexer = StringReversedIndexer();
    auto rootsplit_indices = indexer_representation.first;
    auto pcss_index_vector = indexer_representation.second;
    StringSet rootsplit_string_set;
    for (auto index : rootsplit_indices) {
      rootsplit_string_set.insert(reversed_indexer[index]);
    }
    StringSetVector pcss_string_sets;
    for (const auto &pcss_indices : pcss_index_vector) {
      StringSet pcss_string_set;
      for (auto index : pcss_indices) {
        pcss_string_set.insert(reversed_indexer[index]);
      }
      pcss_string_sets.push_back(std::move(pcss_string_set));
    }
    return {rootsplit_string_set, pcss_string_sets};
  }

  // ** I/O

  std::tuple<StringUInt32Map, StringUInt32PairMap> GetIndexers() {
    auto indexer_str = StringUInt32MapOf(indexer_);
    StringUInt32PairMap parent_to_range_str;
    for (const auto &iter : parent_to_range_) {
      assert(parent_to_range_str.insert({iter.first.ToString(), iter.second})
                 .second);
    }
    assert(parent_to_range_str.insert({"rootsplit", {0, rootsplit_index_end_}})
               .second);
    return std::tie(indexer_str, parent_to_range_str);
  }

  // This function is really just for testing-- it recomputes from scratch.
  std::pair<StringUInt32Map, StringPCSSMap> SplitCounters() {
    auto counter = tree_collection_.TopologyCounter();
    return {StringUInt32MapOf(RootsplitCounterOf(counter)),
            StringPCSSMapOf(PCSSCounterOf(counter))};
  }

  void ReadNewickFile(std::string fname) {
    Driver driver;
    tree_collection_ = driver.ParseNewickFile(fname);
  }

  void ReadNexusFile(std::string fname) {
    Driver driver;
    tree_collection_ = driver.ParseNexusFile(fname);
  }

  void ReadFastaFile(std::string fname) { alignment_.ReadFasta(fname); }

  // ** Phylogenetic likelihood

  void CheckDataLoaded() {
    if (alignment_.SequenceCount() == 0) {
      std::cerr << "Load an alignment into your SBNInstance on which you wish "
                   "to calculate phylogenetic likelihoods.\n";
      abort();
    }
    if (TreeCount() == 0) {
      std::cerr << "Load some trees into your SBNInstance on which you wish to "
                   "calculate phylogenetic likelihoods.\n";
      abort();
    }
  }

  void CheckBeagleDimensions() {
    CheckDataLoaded();
    if (beagle_instances_.size() == 0) {
      std::cerr << "Call MakeBeagleInstances to make some instances for "
                   "likelihood computation.\n";
      abort();
    }
    if (alignment_.SequenceCount() != beagle_leaf_count_ ||
        alignment_.Length() != beagle_site_count_) {
      std::cerr << "Alignment dimensions for current BEAGLE instances do not "
                   "match current alignment. Call MakeBeagleInstances again.\n";
      abort();
    }
  }

  void MakeBeagleInstances(int instance_count) {
    // Start by clearing out any existing instances.
    FinalizeBeagleInstances();
    CheckDataLoaded();
    beagle_leaf_count_ = alignment_.SequenceCount();
    beagle_site_count_ = alignment_.Length();
    SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
    for (auto i = 0; i < instance_count; i++) {
      auto beagle_instance = beagle::CreateInstance(site_pattern);
      beagle::SetJCModel(beagle_instance);
      beagle_instances_.push_back(beagle_instance);
      beagle::PrepareBeagleInstance(beagle_instance, tree_collection_,
                                    site_pattern);
    }
  }

  std::vector<double> LogLikelihoods() {
    CheckBeagleDimensions();
    return beagle::LogLikelihoods(beagle_instances_, tree_collection_,
                                  rescaling_);
  }

  std::vector<std::pair<double, std::vector<double>>> BranchGradients() {
    return beagle::BranchGradients(beagle_instances_, tree_collection_,
                                   rescaling_);
  }
};

// Here we initialize our static random number generator.
std::random_device SBNInstance::random_device_;
std::mt19937 SBNInstance::random_generator_(SBNInstance::random_device_());

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
  // See https://github.com/matsengrp/libsbn/issues/74 to understand this test.
  auto indexer_test_topology_1 =
      // (2,(1,3)5,(0,4)6)7
      Node::OfParentIndexVector({6, 5, 7, 5, 6, 7, 7});
  std::cout << indexer_test_topology_1->Newick(
      [](const Node *node) { return std::to_string(node->Index()); });
  IndexerRepresentation correct_representation_1(
      {{8, 0, 3, 5, 10, 4, 11},  // The rootsplit indices.
       {{54, 31, 15},            // The PCSS indices.
        {21, 13, 22},
        {74, 42, 61},
        {21, 12, 63},
        {27, 31, 57},
        {21, 75, 76},
        {70, 73, 31}}});
  std::cout << inst.StringIndexerRepresentationOf(correct_representation_1)
            << std::endl;
  CHECK_EQ(IndexerRepresentationOf(inst.indexer_, indexer_test_topology_1),
           correct_representation_1);
  auto indexer_test_topology_2 =
      // (((0,1)5,2)6,3,4)7;
      Node::OfParentIndexVector({5, 5, 6, 7, 7, 6, 7});
  IndexerRepresentation correct_representation_2(
      {{8, 0, 3, 5, 10, 6, 1},  // The rootsplit indices.
       {{50, 36, 16},           // The PCSS indices.
        {50, 49, 23},
        {39, 34, 44},
        {29, 65, 72},
        {59, 41, 72},
        {50, 47, 51},
        {40, 30, 72}}});
  std::cout << inst.StringIndexerRepresentationOf(correct_representation_2)
            << std::endl;
  CHECK_EQ(IndexerRepresentationOf(inst.indexer_, indexer_test_topology_2),
           correct_representation_2);
  inst.SampleTrees(2);
  inst.GetIndexerRepresentations();

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
    std::cout << gradients[i].first << " " << pybeagle_likelihoods[i]
              << std::endl;
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
