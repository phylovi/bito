// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIBSBN_HPP_
#define SRC_LIBSBN_HPP_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <cmath>
#include <queue>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "alignment.hpp"
#include "beagle.hpp"
#include "build.hpp"
#include "driver.hpp"
#include "task_processor.hpp"
#include "tree.hpp"

namespace py = pybind11;

typedef std::unordered_map<std::string, float> StringFloatMap;
typedef std::unordered_map<std::string, uint32_t> StringUInt32Map;

template <typename T>
StringUInt32Map StringUInt32MapOf(T m) {
  StringUInt32Map m_str;
  for (const auto &iter : m) {
    m_str[iter.first.ToString()] = iter.second;
  }
  return m_str;
}

template <typename T>
BitsetUInt32Map BitsetIndexerOf(T supports) {
  BitsetUInt32Map indexer;
  uint32_t index = 0;
  for (const auto &iter : supports) {
    indexer[iter.first] = index;
    index++;
  }
  return indexer;
}

struct SBNInstance {
  std::string name_;
  TreeCollection::TreeCollectionPtr tree_collection_;
  Alignment alignment_;
  CharIntMap symbol_table_;
  std::vector<beagle::BeagleInstance> beagle_instances_;
  std::vector<double> sbn_probs_;
  BitsetUInt32Map rootsplit_indexer_;
  BitsetUInt32Map pcss_indexer_;

  // ** Initialization, destruction, and status
  explicit SBNInstance(const std::string &name)
      : name_(name), symbol_table_(beagle::GetSymbolTable()) {}

  ~SBNInstance() { FinalizeBeagleInstances(); }

  // Finalize means to release memory.
  void FinalizeBeagleInstances() {
    for (const auto &beagle_instance : beagle_instances_) {
      assert(beagleFinalizeInstance(beagle_instance) == 0);
    }
    beagle_instances_.clear();
  }

  size_t TreeCount() const { return tree_collection_->TreeCount(); }
  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (tree_collection_) {
      std::cout << TreeCount() << " unique tree topologies loaded on "
                << tree_collection_->TaxonCount() << " leaves.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
    std::cout << alignment_.Data().size() << " sequences loaded.\n";
  }

  // ** Building SBN-related items

  void ProcessLoadedTrees() {
    auto counter = tree_collection_->TopologyCounter();
    rootsplit_indexer_ = BitsetIndexerOf(RootsplitCounterOf(counter));
    pcss_indexer_ = BitsetIndexerOf(PCSSCounterOf(counter));
    sbn_probs_ = std::vector<double>(rootsplit_indexer_.size());
  }

  // TODO(erick) replace with something interesting.
  double SBNTotalProb() {
    double total = 0;
    for (const auto &prob : sbn_probs_) {
      total += prob;
    }
    return total;
  }

  // ** I/O

  StringUInt32Map GetRootsplitIndexer() {
    return StringUInt32MapOf(rootsplit_indexer_);
  }

  StringUInt32Map GetPCSSIndexer() { return StringUInt32MapOf(pcss_indexer_); }

  // This function is really just for testing-- it recomputes from scratch.
  std::pair<StringUInt32Map, StringUInt32Map> SplitCounters() {
    auto counter = tree_collection_->TopologyCounter();
    return {StringUInt32MapOf(RootsplitCounterOf(counter)),
            StringUInt32MapOf(PCSSCounterOf(counter))};
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

  void MakeBeagleInstances(int instance_count) {
    // Start by clearing out any existing instances.
    FinalizeBeagleInstances();
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
    for (auto i = 0; i < instance_count; i++) {
      auto beagle_instance = beagle::CreateInstance(alignment_);
      beagle::SetJCModel(beagle_instance);
      beagle_instances_.push_back(beagle_instance);
      beagle::PrepareBeagleInstance(beagle_instance, tree_collection_,
                                    alignment_, symbol_table_);
    }
  }

  std::vector<double> TreeLogLikelihoods() {
    if (beagle_instances_.size() == 0) {
      std::cerr << "Please add some BEAGLE instances that can be used for "
                   "computation.\n";
      abort();
    }
    std::vector<double> results(tree_collection_->TreeCount());
    std::queue<beagle::BeagleInstance> instance_queue;
    for (auto instance : beagle_instances_) {
      instance_queue.push(instance);
    }
    std::queue<size_t> tree_number_queue;
    for (size_t i = 0; i < tree_collection_->TreeCount(); i++) {
      tree_number_queue.push(i);
    }
    TaskProcessor<beagle::BeagleInstance, size_t> task_processor(
        instance_queue, tree_number_queue,
        [&results, &tree_collection = tree_collection_ ](
            beagle::BeagleInstance beagle_instance, size_t tree_number) {
          results[tree_number] = beagle::TreeLogLikelihood(
              tree_collection->GetTree(tree_number), beagle_instance);
        });
    return results;
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("libsbn") {
  SBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon.nwk");
  // Reading one file after another checks that we've cleared out state.
  inst.ReadNewickFile("data/hello.nwk");
  inst.ReadFastaFile("data/hello.fasta");
  inst.MakeBeagleInstances(2);
  for (auto ll : inst.TreeLogLikelihoods()) {
    CHECK_LT(abs(ll - -84.852358), 0.000001);
  }
  inst.ReadNexusFile("data/DS1.subsampled_10.t");
  inst.ReadFastaFile("data/DS1.fasta");
  inst.MakeBeagleInstances(2);
  auto likelihoods = inst.TreeLogLikelihoods();
  std::vector<double> pybeagle_likelihoods(
      {-14582.995273982739, -6911.294207416366, -6916.880235529542,
       -6904.016888831189, -6915.055570693576, -6915.50496696512,
       -6910.958836661867, -6909.02639968063, -6912.967861935749,
       -6910.7871105783515});
  for (size_t i = 0; i < likelihoods.size(); i++) {
    CHECK_LT(abs(likelihoods[i] - pybeagle_likelihoods[i]), 0.00011);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_LIBSBN_HPP_
