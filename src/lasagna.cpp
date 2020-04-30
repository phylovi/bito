
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "gp_instance.hpp"

TEST_CASE("GPInstance") {
  // ((mars:0.1,saturn:0.1):0.2,jupiter:0.1):0.1;
  GPOperationVector operations{
      GPOperations::SetToStationaryDistribution{3},
      // Evolve mars rootward.
      GPOperations::EvolveRootward{4, 0, 0},
      // Evolve saturn rootward.
      GPOperations::EvolveRootward{5, 1, 1},
      // Get common ancestor of mars and saturn.
      GPOperations::Multiply{6, 4, 5},
      // Evolve common ancestor rootward.
      GPOperations::EvolveRootward{7, 6, 2},
      // Evolve jupiter rootward.
      GPOperations::EvolveRootward{8, 2, 3},
      // Get root.
      GPOperations::Multiply{9, 7, 8},
      // Calculate likelihood.
      GPOperations::Likelihood{0, 3, 9},
  };

  GPInstance inst;
  inst.ReadFastaFile("data/hello.fasta");
  inst.ReadNewickFile("data/hello_rooted.nwk");
  inst.MakeEngine();
  auto engine = inst.GetEngine();
  EigenVectorXd branch_lengths(4);
  branch_lengths << 0.1, 0.1, 0.2, 0.1;
  engine->SetBranchLengths(branch_lengths);
  engine->ProcessOperations(operations);
  CHECK_LT(fabs(engine->GetLikelihood(0) - -84.852358), 1e-6);
}
