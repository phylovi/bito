
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "gp_instance.hpp"

GPInstance MakeHelloGPInstance() {
  // ((mars:0.1,saturn:0.1):0.2,jupiter:0.1):0.1;
  GPInstance inst;
  inst.ReadFastaFile("data/hello.fasta");
  inst.ReadNewickFile("data/hello_rooted.nwk");
  inst.MakeEngine();
  EigenVectorXd branch_lengths(4);
  branch_lengths << 0.1, 0.1, 0.2, 0.1;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  return inst;
}

GPOperationVector rootward_likelihood_calculation{
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

TEST_CASE("GPInstance: rootward likelihood calculation") {
  auto inst = MakeHelloGPInstance();
  auto engine = inst.GetEngine();
  engine->ProcessOperations(rootward_likelihood_calculation);
  CHECK_LT(fabs(engine->GetLogLikelihood(0) - -84.852358), 1e-6);
  }
