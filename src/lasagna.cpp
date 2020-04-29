#include <vector>
#include "gp_operation.hpp"

int main() {
  // auto z = GPOperations::Zero{5};
  auto z = GPOperations::WeightedSumAccumulate{5, 4, 5};
  std::cout << z << std::endl;
}
