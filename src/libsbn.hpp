// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.


#ifndef SRC_LIBSBN_HPP_
#define SRC_LIBSBN_HPP_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <string>
#include <unordered_map>
#include "alignment.hpp"
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

  explicit SBNInstance(const std::string &name) : name_(name) {}

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
