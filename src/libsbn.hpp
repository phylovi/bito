#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <unordered_map>
#include "build.hpp"
#include "driver.hpp"
#include "tree.hpp"

namespace py = pybind11;

typedef std::unordered_map<std::string, float> StringFloatMap;

StringFloatMap StringFloatMapOf(BitsetFloatMap m) {
  StringFloatMap m_str;
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    m_str[iter->first.ToString()] = iter->second;
  }
  return m_str;
}


struct SBNInstance {
  std::string name_;
  Driver driver_;
  Node::NodePtrVecPtr trees_;
  std::unordered_map<int, int> indexer_;

  SBNInstance(const std::string &name) : name_(name) {}

  unsigned long TreeCount() { return trees_->size(); }
  void ParseFile(std::string fname) { trees_ = driver_.ParseFile(fname); }

  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (trees_) {
      std::cout << trees_->size() << " trees loaded.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
  }

  std::unordered_map<uint64_t, std::string> g() {
    assert(trees_->size() > 0);
    TagBitsetMap m = TagBitsetMapOf(trees_->at(0));
    std::unordered_map<uint64_t, std::string> m_str;
    m = TagBitsetMapOf(trees_->at(0));
    for (auto iter = m.begin(); iter != m.end(); ++iter) {
      m_str[iter->first] = iter->second.ToString();
    }
    return m_str;
  }

  StringFloatMap Rootsplits() {
    assert(trees_->size() > 0);
    return StringFloatMapOf(RootsplitFrequencyOf(trees_->at(0)));
  }

  static void f(py::array_t<double> array) {
    py::buffer_info buf = array.request();
    std::cout << "You passed a " << buf.ndim << " dim array" << std::endl;
    double *ptr = (double *)buf.ptr;
    for (auto idx = 0; idx < buf.shape[0]; idx++) {
      std::cout << ptr[idx] << std::endl;
    }
    };
};
