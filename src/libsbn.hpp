#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <unordered_map>
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
    m_str[iter->first.ToPCSSString()] = iter->second;
  }
  return m_str;
}

struct SBNInstance {
  std::string name_;
  Driver driver_;
  Node::NodePtrCounterPtr trees_;
  std::unordered_map<int, int> indexer_;

  SBNInstance(const std::string &name) : name_(name) {}

  unsigned long TreeCount() { return trees_->size(); }
  void ParseFile(std::string fname) { trees_ = driver_.ParseFile(fname); }

  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (trees_) {
      std::cout << trees_->size() << " unique tree topologies loaded.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
  }

  StringUInt32Map Supports() { return StringUInt32MapOf(SupportsOf(trees_)); }

  static void f(py::array_t<double> array) {
    py::buffer_info buf = array.request();
    std::cout << "You passed a " << buf.ndim << " dim array" << std::endl;
    double *ptr = (double *)buf.ptr;
    for (auto idx = 0; idx < buf.shape[0]; idx++) {
      std::cout << ptr[idx] << std::endl;
    }
    };
};
