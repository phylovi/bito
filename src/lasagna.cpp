
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

// Here "leaf" is the leaf-side PLV, and "root" is the root-side PLV.
enum HelloPLV {
  mars_leaf = 0,
  saturn_leaf,
  jupiter_leaf,
  mars_root,
  saturn_root,
  jupiter_root,
  ancestor_leaf,
  ancestor_root,
  root,
  stationary,
};

enum HelloBranchLength { mars, saturn, ancestor, jupiter };

GPOperationVector rootward_likelihood_calculation{
    GPOperations::SetToStationaryDistribution{HelloPLV::stationary},
    // Evolve mars rootward.
    GPOperations::EvolveRootward{HelloPLV::mars_root, HelloPLV::mars_leaf,
                                 HelloBranchLength::mars},
    // Evolve saturn rootward.
    GPOperations::EvolveRootward{HelloPLV::saturn_root, HelloPLV::saturn_leaf,
                                 HelloBranchLength::saturn},
    // Get common ancestor of mars and saturn.
    GPOperations::Multiply{HelloPLV::ancestor_leaf, HelloPLV::mars_root,
                           HelloPLV::saturn_root},
    // Evolve common ancestor rootward.
    GPOperations::EvolveRootward{HelloPLV::ancestor_root, HelloPLV::ancestor_leaf,
                                 HelloBranchLength::ancestor},
    // Evolve jupiter rootward.
    GPOperations::EvolveRootward{HelloPLV::jupiter_root, HelloPLV::jupiter_leaf,
                                 HelloBranchLength::jupiter},
    // Get root.
    GPOperations::Multiply{HelloPLV::root, HelloPLV::ancestor_root,
                           HelloPLV::jupiter_root},
    // Calculate likelihood.
    GPOperations::Likelihood{0, HelloPLV::stationary, HelloPLV::root},
};

TEST_CASE("GPInstance: rootward likelihood calculation") {
  auto inst = MakeHelloGPInstance();
  auto engine = inst.GetEngine();
  engine->ProcessOperations(rootward_likelihood_calculation);
  CHECK_LT(fabs(engine->GetLogLikelihoods()(0) - -84.852358), 1e-6);
}

TEST_CASE("GPInstance: leafward likelihood calculation") {
  auto inst = MakeHelloGPInstance();
  auto engine = inst.GetEngine();
  engine->ProcessOperations(rootward_likelihood_calculation);
  // engine->ProcessOperations(leafward_likelihood_calculation);
  CHECK_LT(fabs(engine->GetLogLikelihoods()(0) - -84.852358), 1e-6);
}
