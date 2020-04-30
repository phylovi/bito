#include "gp_operation.hpp"

int main() {
  GPOperationVector operations{GPOperations::Zero{5},
                               GPOperations::UpdateSBNProbabilities{2, 5}};

  std::cout << operations << std::endl;
}
