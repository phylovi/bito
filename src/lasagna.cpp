
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "gp_instance.hpp"

using namespace GPOperations;

enum HelloGPCSP { jupiter, mars, saturn, ancestor, root };

// According to HelloGPCSP, this makes
// (jupiter:0.113,(mars:0.15,saturn:0.1):0.22):0.;
GPInstance MakeHelloGPInstance() {
  GPInstance inst;
  inst.ReadFastaFile("data/hello.fasta");
  inst.ReadNewickFile("data/hello_rooted.nwk");
  inst.MakeEngine();
  EigenVectorXd branch_lengths(4);
  // Order set by HelloGPCSP.
  branch_lengths << 0.113, 0.15, 0.1, 0.22;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  return inst;
}

TEST_CASE("GPInstance: straightforward classical likelihood calculation") {
  // Here "_leaf" is the leaf-side PLV, and "_root" is the root-side PLV.
  enum PLV {
    jupiter_leaf,
    mars_leaf,
    saturn_leaf,
    jupiter_root,
    mars_root,
    saturn_root,
    ancestor_leaf,
    ancestor_root,
    root_leaf,  // The leaf-side PLV for the root node.
    stationary,
  };

  GPOperationVector rootward_likelihood_calculation{
      SetToStationaryDistribution{PLV::stationary},
      // Evolve mars rootward.
      EvolveRootward{PLV::mars_root, PLV::mars_leaf, HelloGPCSP::mars},
      // Evolve saturn rootward.
      EvolveRootward{PLV::saturn_root, PLV::saturn_leaf, HelloGPCSP::saturn},
      // Get common ancestor of mars and saturn.
      Multiply{PLV::ancestor_leaf, PLV::mars_root, PLV::saturn_root},
      // Evolve common ancestor rootward.
      EvolveRootward{PLV::ancestor_root, PLV::ancestor_leaf, HelloGPCSP::ancestor},
      // Evolve jupiter rootward.
      EvolveRootward{PLV::jupiter_root, PLV::jupiter_leaf, HelloGPCSP::jupiter},
      // Get root.
      Multiply{PLV::root_leaf, PLV::ancestor_root, PLV::jupiter_root},
      // Calculate likelihood at the root.
      Likelihood{HelloGPCSP::root, PLV::stationary, PLV::root_leaf},
  };

  auto inst = MakeHelloGPInstance();
  auto engine = inst.GetEngine();
  engine->ProcessOperations(rootward_likelihood_calculation);
  CHECK_LT(fabs(engine->GetLogLikelihoods()(HelloGPCSP::root) - -84.77961943), 1e-6);
}

TEST_CASE("GPInstance: subsplit traversal as written") {
  // Here t is the rootsplit, and s is the subsplit with mars and jupiter.
  enum PLV { p_ };

  /*
    GPOperationVector leafward_likelihood_calculation{
        // Get the PLV heading down towards the ancestor.
        Multiply{HelloPLV::ancestor_root, HelloPLV::stationary, HelloPLV::jupiter_root},
        // Evolve common ancestor leafward.
        EvolveLeafward{HelloPLV::ancestor_leaf, HelloPLV::ancestor_root,
                       HelloGPCSP::ancestor},
        // Calculate likelihood at common ancestor.
        Likelihood{HelloGPCSP::ancestor, HelloPLV::ancestor_leaf, HelloPLV::root_leaf},
        // Get common ancestor of mars and saturn.
        Multiply{HelloPLV::ancestor_leaf, HelloPLV::mars_root, HelloPLV::saturn_root},
        // Evolve common ancestor rootward.
        // Evolve jupiter rootward.
        EvolveRootward{HelloPLV::jupiter_root, HelloPLV::jupiter_leaf,
                       HelloGPCSP::jupiter},
        // Get root.
        Multiply{HelloPLV::root_leaf, HelloPLV::ancestor_root, HelloPLV::jupiter_root},
        // Calculate likelihood.
        Likelihood{HelloGPCSP::root, HelloPLV::stationary, HelloPLV::root_leaf},
    };

    auto inst = MakeHelloGPInstance();
    auto engine = inst.GetEngine();
    engine->ProcessOperations(rootward_likelihood_calculation);
    // engine->ProcessOperations(leafward_likelihood_calculation);
    // CHECK_LT(fabs(engine->GetLogLikelihoods()(0) - -84.852358), 1e-6);
  */
}
