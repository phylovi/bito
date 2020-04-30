
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "gp_instance.hpp"

TEST_CASE("GPInstance") {
  GPOperationVector operations{GPOperations::Zero{5},
                               GPOperations::UpdateSBNProbabilities{2, 5}};

  std::cout << operations << std::endl;

  GPInstance inst;

  inst.ReadFastaFile("data/hello.fasta");
  inst.ReadNewickFile("data/hello_rooted.nwk");
  inst.MakeEngine();
}
