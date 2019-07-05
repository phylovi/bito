// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIBSBN_HPP_
#define SRC_LIBSBN_HPP_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <string>
#include <unordered_map>
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
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    m_str[iter->first.ToString()] = static_cast<float>(iter->second);
  }
  return m_str;
}

StringUInt32Map StringUInt32MapOf(BitsetUInt32Map m) {
  StringUInt32Map m_str;
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    m_str[iter->first.ToString()] = iter->second;
  }
  return m_str;
}

struct SBNInstance {
  std::string name_;
  Driver driver_;
  TreeCollection::TreeCollectionPtr tree_collection_;
  Alignment alignment_;
  CharIntMap symbol_table_ = GetSymbolTable();
  int beagle_instance_ = -1;

  explicit SBNInstance(const std::string &name) : name_(name) {
  }

  ~SBNInstance() {}

  size_t TreeCount() { return tree_collection_->TreeCount(); }

  void ParseFile(std::string fname) {
    tree_collection_ = driver_.ParseFile(fname);
  }

  void ReadFasta(std::string fname) { alignment_.ReadFasta(fname); }

  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (tree_collection_) {
      std::cout << TreeCount() << " unique tree topologies loaded.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
    std::cout << alignment_.Data().size() << " sequences loaded.\n";
  }

  StringUInt32Map RootsplitSupport() {
    return StringUInt32MapOf(RootsplitSupportOf(tree_collection_->Trees()));
  }
  StringUInt32Map SubsplitSupport() {
    return StringUInt32MapOf(SubsplitSupportOf(tree_collection_->Trees()));
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
    int tip_count = tree_collection_->FirstTree()->LeafCount();
    if (tip_count != int(alignment_.SequenceCount())) {
      std::cerr << "The number of tree tips doesn't match the alignment "
                   "sequence count!\n";
      abort();
    }
    // Number of partial buffers to create (input) -- internal node count
    int partials_buffer_count = tip_count - 1;
    // Number of compact state representation buffers to create -- for use with
    // setTipStates (input) */
    int compact_buffer_count = tip_count;
    // DNA assumption here.
    int state_count = 4;
    // Number of site patterns to be handled by the instance (input) -- not
    // compressed in this case
    int pattern_count = int(alignment_.Length());
    // Number of eigen-decomposition buffers to allocate (input)
    int eigen_buffer_count = 1;
    // Number of transition matrix buffers (input) -- one per edge
    int matrix_buffer_count = 2 * tip_count - 1;
    // Number of rate categories
    int category_count = 1;
    // Number of scaling buffers -- can be zero if scaling is not needed
    int scale_buffer_count = 0;
    // List of potential resources on which this instance is allowed (input,
    // NULL implies no restriction
    int *allowed_resources = NULL;
    // Length of resourceList list (input) -- not needed to use the default
    // hardware config
    int resource_count = 0;
    // Bit-flags indicating preferred implementation charactertistics, see
    // BeagleFlags (input)
    long preference_flags = 0;
    // Bit-flags indicating required implementation characteristics, see
    // BeagleFlags (input)
    int requirement_flags = 0;

    BeagleInstanceDetails *return_info = new BeagleInstanceDetails();

    beagle_instance_ = beagleCreateInstance(
        tip_count, partials_buffer_count, compact_buffer_count, state_count,
        pattern_count, eigen_buffer_count, matrix_buffer_count, category_count,
        scale_buffer_count, allowed_resources, resource_count, preference_flags,
        requirement_flags, return_info);

    delete return_info;
  }

  static void f(py::array_t<double> array) {
    py::buffer_info buf = array.request();
    std::cout << "You passed a " << buf.ndim << " dim array" << std::endl;
    double *ptr = (double *)buf.ptr;
    for (auto idx = 0; idx < buf.shape[0]; idx++) {
      std::cout << ptr[idx] << std::endl;
    }
  }
};

#endif  // SRC_LIBSBN_HPP_
