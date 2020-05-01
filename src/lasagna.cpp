
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "gp_instance.hpp"

using namespace GPOperations;

// The "ancestor" node is the common ancestor of mars and saturn.
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
  // Here t is the rootsplit (jupiter, other), and s is the osubsplit (mars, jupiter).
  enum PLV {
    p_jupiter,     // fake subsplit leading to jupiter
    p_mars,        // fake subsplit leading to mars
    p_saturn,      // fake subsplit leading to mars
    phat_stilde,   // root side of mars edge
    phat_s,        // root side of saturn edge
    p_s,           // p at ancestor
    phat_ttilde,   // root side of jupiter edge
    phat_t,        // root side of ancestor edge
    p_t,           // PLV at root from leaves
    rhat_t,        // stationary distribution coming from root
    r_ttilde,      // PLV pointing towards jupiter from root
    r_t,           // PLV pointing towards ancestor from root
    rhat_s,        // leaf side of ancestor edge
    r_stilde,      // PLV pointing towards mars from ancestor
    r_s,           // PLV pointing towards saturn from ancestor
    rhat_jupiter,  // Evolved message coming towards jupiter
    rhat_mars,     // Evolved message coming towards mars
    rhat_saturn,   // Evolved message coming towards saturn
  };

  GPOperationVector rootward_likelihood_calculation{
      // Evolve jupiter rootward.
      EvolveRootward{PLV::phat_ttilde, PLV::p_jupiter, HelloGPCSP::jupiter},
      // Evolve mars rootward.
      EvolveRootward{PLV::phat_stilde, PLV::p_mars, HelloGPCSP::mars},
      // Evolve saturn rootward.
      EvolveRootward{PLV::phat_s, PLV::p_saturn, HelloGPCSP::saturn},
      // Get common ancestor of mars and saturn.
      Multiply{PLV::p_s, PLV::phat_stilde, PLV::phat_s},
      // Evolve ancestor rootward.
      EvolveRootward{PLV::phat_t, PLV::p_s, HelloGPCSP::ancestor},
      // Get root PLV coming from leaves
      Multiply{PLV::p_t, PLV::phat_ttilde, PLV::phat_t},
      // Set a stationary distribution coming from "beyond" the root.
      SetToStationaryDistribution{PLV::rhat_t},
      // Calculate likelihood at root.
      Likelihood{HelloGPCSP::root, PLV::rhat_t, PLV::p_t},
      // Get message going towards jupiter.
      Multiply{PLV::r_ttilde, PLV::rhat_t, PLV::phat_t},
      // Evolve the message towards jupiter.
      EvolveLeafward{PLV::rhat_jupiter, PLV::r_ttilde, HelloGPCSP::jupiter},
      // Calculate likelihood at jupiter.
      Likelihood{HelloGPCSP::jupiter, PLV::rhat_jupiter, PLV::p_jupiter},
      // Get message going towards ancestor.
      Multiply{PLV::r_t, PLV::phat_ttilde, PLV::rhat_t},
      // Evolve the message towards ancestor.
      EvolveLeafward{PLV::rhat_s, PLV::r_t, HelloGPCSP::ancestor},
      // Calculate likelihood at common ancestor.
      Likelihood{HelloGPCSP::ancestor, PLV::rhat_s, PLV::p_s},
      // Get message going towards mars.
      Multiply{PLV::r_stilde, PLV::rhat_s, PLV::phat_s},
      // Evolve the message towards mars.
      EvolveLeafward{PLV::rhat_mars, PLV::r_stilde, HelloGPCSP::mars},
      // Calculate likelihood at mars.
      Likelihood{HelloGPCSP::mars, PLV::rhat_mars, PLV::p_mars},
      // Get message going towards saturn.
      Multiply{PLV::r_s, PLV::rhat_s, PLV::phat_stilde},
      // Evolve the message towards saturn.
      EvolveLeafward{PLV::rhat_saturn, PLV::r_s, HelloGPCSP::saturn},
      // Calculate likelihood at saturn.
      Likelihood{HelloGPCSP::saturn, PLV::rhat_saturn, PLV::p_saturn},

  };

  auto inst = MakeHelloGPInstance();
  auto engine = inst.GetEngine();
  engine->ProcessOperations(rootward_likelihood_calculation);
  CHECK_LT(fabs(engine->GetLogLikelihoods()(HelloGPCSP::root) - -84.77961943), 1e-6);
  CHECK_LT(fabs(engine->GetLogLikelihoods()(HelloGPCSP::mars) - -84.77961943), 1e-6);
  CHECK_LT(fabs(engine->GetLogLikelihoods()(HelloGPCSP::jupiter) - -84.77961943), 1e-6);
  CHECK_LT(fabs(engine->GetLogLikelihoods()(HelloGPCSP::saturn) - -84.77961943), 1e-6);
  CHECK_LT(fabs(engine->GetLogLikelihoods()(HelloGPCSP::ancestor) - -84.77961943),
           1e-6);
}
