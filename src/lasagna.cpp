#include <vector>
#include "gp_operation.hpp"

int main() {
  // GPOperationVector z{GPOperations::Zero{5},
  //                     GPOperations::WeightedSumAccumulate{5, 4, 5}};
  GPOperation z = GPOperations::Zero{5};
  std::cout << z << std::endl;
}
