
#include <unordered_map>
#include "driver.hpp"
#include "tree.hpp"

struct SBNInstance {
  std::string name_;
  Driver driver_;
  Node::NodePtrVecPtr trees_;
  std::unordered_map<int, int> indexer_;

  SBNInstance(const std::string &name) : name_(name) {}

  unsigned long TreeCount() { return trees_->size(); }
  void ParseFile(std::string fname) { trees_ = driver_.ParseFile(fname); }

  void InitIndexer() { indexer_[4] = 2; }

  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (trees_) {
      std::cout << trees_->size() << " trees loaded.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
  }
};
