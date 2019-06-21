
#include "driver.hpp"
#include "tree.hpp"

class SBNInstance {
 public:
  Driver driver_;
  Node::NodePtrVecPtr trees_;

  void ParseFile(std::string fname) { trees_ = driver_.ParseFile(fname); }

  void PrintStatus() {
    if (trees_) {
      std::cout << trees_->size() << " trees loaded.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
  }
};
