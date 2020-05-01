
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "gp_instance.hpp"

using namespace GPOperations;

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

enum HelloGPCSP { mars, saturn, ancestor, jupiter, root };

// Here "_leaf" is the leaf-side PLV, and "_root" is the root-side PLV.
enum HelloPLV {
  mars_leaf = 0,
  saturn_leaf,
  jupiter_leaf,
  mars_root,
  saturn_root,
  jupiter_root,
  ancestor_leaf,
  ancestor_root,
  root_leaf,  // The leaf-side PLV for the root node.
  stationary,
};

GPOperationVector rootward_likelihood_calculation{
    SetToStationaryDistribution{stationary},
    // Evolve mars rootward.
    EvolveRootward{HelloPLV::mars_root, HelloPLV::mars_leaf, HelloGPCSP::mars},
    // Evolve saturn rootward.
    EvolveRootward{HelloPLV::saturn_root, HelloPLV::saturn_leaf, HelloGPCSP::saturn},
    // Get common ancestor of mars and saturn.
    Multiply{HelloPLV::ancestor_leaf, HelloPLV::mars_root, HelloPLV::saturn_root},
    // Evolve common ancestor rootward.
    EvolveRootward{HelloPLV::ancestor_root, HelloPLV::ancestor_leaf,
                   HelloGPCSP::ancestor},
    // Evolve jupiter rootward.
    EvolveRootward{HelloPLV::jupiter_root, HelloPLV::jupiter_leaf, HelloGPCSP::jupiter},
    // Get root.
    Multiply{HelloPLV::root_leaf, HelloPLV::ancestor_root, HelloPLV::jupiter_root},
    // Calculate likelihood.
    Likelihood{0, HelloPLV::stationary, HelloPLV::root_leaf},
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
