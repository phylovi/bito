// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIBSBN_HPP_
#define SRC_LIBSBN_HPP_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <cmath>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "alignment.hpp"
#include "beagle.hpp"
#include "build.hpp"
#include "driver.hpp"
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
  beagle::BeagleInstance beagle_instance_;

  explicit SBNInstance(const std::string &name)
      : name_(name),
        symbol_table_(beagle::GetSymbolTable()),
        beagle_instance_(-1) {}

  ~SBNInstance() { assert(beagleFinalizeInstance(beagle_instance_) == 0); }

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

  void BeagleCreate() {
    if (beagle_instance_ != -1) {
      std::cerr << "BEAGLE Instance already exists and needs to be freed.\n";
      abort();
    }
    if (TreeCount() == 0) {
      std::cerr << "Load some trees into your SBNInstance on which you wish to "
                   "calculate phylogenetic likelihoods.\n";
      abort();
    }
    if (alignment_.SequenceCount() == 0) {
      std::cerr << "Load an alignment into your SBNInstance on which you wish "
                   "to calculate phylogenetic likelihoods.\n";
      abort();
    }
    beagle_instance_ = beagle::CreateInstance(tree_collection_, alignment_);
  }

  void PrepareBeagleInstance() {
    assert(beagle_instance_ != -1);
    beagle::PrepareBeagleInstance(beagle_instance_, tree_collection_,
                                  alignment_, symbol_table_);
  }

  void SetJCModel() {
    assert(beagle_instance_ != -1);
    beagle::SetJCModel(beagle_instance_);
  }

  std::vector<double> TreeLogLikelihoods() {
    assert(beagle_instance_ != -1);
    return beagle::TreeLogLikelihoods(beagle_instance_, tree_collection_);
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
  inst.BeagleCreate();
  inst.PrepareBeagleInstance();
  inst.SetJCModel();
  for (auto ll : inst.TreeLogLikelihoods()) {
    CHECK_LT(abs(ll - -84.852358), 0.000001);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_LIBSBN_HPP_
