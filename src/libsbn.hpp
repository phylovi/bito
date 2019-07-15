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

StringFloatMap StringFloatMapOf(BitsetUInt32Map m) {
  StringFloatMap m_str;
  for (const auto &iter : m) {
    m_str[iter.first.ToString()] = static_cast<float>(iter.second);
  }
  return m_str;
}

StringUInt32Map StringUInt32MapOf(BitsetUInt32Map m) {
  StringUInt32Map m_str;
  for (const auto &iter : m) {
    m_str[iter.first.ToString()] = iter.second;
  }
  return m_str;
}

struct SBNInstance {
  std::string name_;
  TreeCollection::TreeCollectionPtr tree_collection_;
  Alignment alignment_;
  CharIntMap symbol_table_;
  std::vector<beagle::BeagleInstance> beagle_instances_;

  explicit SBNInstance(const std::string &name)
      : name_(name), symbol_table_(beagle::GetSymbolTable()) {}

  ~SBNInstance() {
    for (const auto &beagle_instance : beagle_instances_) {
      assert(beagleFinalizeInstance(beagle_instance) == 0);
    }
  }

  size_t TreeCount() const { return tree_collection_->TreeCount(); }

  void ReadNewickFile(std::string fname) {
    Driver driver;
    tree_collection_ = driver.ParseNewickFile(fname);
  }

  void ReadFastaFile(std::string fname) { alignment_.ReadFasta(fname); }

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

  std::pair<StringUInt32Map, StringUInt32Map> SplitSupports() {
    auto counter = tree_collection_->TopologyCounter();
    return {StringUInt32MapOf(RootsplitSupportOf(counter)),
            StringUInt32MapOf(SubsplitSupportOf(counter))};
  }

  void AddBeagleInstances(int instance_count) {
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

  static void f(py::array_t<double> array) {
    py::buffer_info buf = array.request();
    std::cout << "You passed a " << buf.ndim << " dim array" << std::endl;
    double *ptr = reinterpret_cast<double *>(buf.ptr);
    for (auto idx = 0; idx < buf.shape[0]; idx++) {
      std::cout << ptr[idx] << std::endl;
    }
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("libsbn") {
  SBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon.nwk");
  inst.PrintStatus();
  inst.ReadNewickFile("data/hello.nwk");
  std::cout << inst.tree_collection_->Newick();
  inst.PrintStatus();
  inst.ReadFastaFile("data/hello.fasta");
  inst.AddBeagleInstances(2);
  for (auto ll : inst.TreeLogLikelihoods()) {
    CHECK_LT(abs(ll - -84.852358), 0.000001);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_LIBSBN_HPP_
