#include "gp_operation.hpp"

int main() {
  GPOperationVector z{GPOperations::Zero{5},
                      GPOperations::UpdateSBNProbabilities{2, 5}};
  std::cout << z << std::endl;
}
