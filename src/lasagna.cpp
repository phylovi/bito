#include "gp_instance.hpp"

int main() {
  GPOperationVector operations{GPOperations::Zero{5},
                               GPOperations::UpdateSBNProbabilities{2, 5}};

  std::cout << operations << std::endl;

  GPInstance inst;

  inst.ReadFastaFile("data/hello.fasta");
  inst.ReadNewickFile("data/hello_rooted.nwk");
  inst.MakeEngine();
  inst.ProcessLoadedTrees();
}
