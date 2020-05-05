
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "gp_instance.hpp"

using namespace GPOperations;

// GPCSP stands for generalized PCSP-- see text.
// Let the "venus" node be the common ancestor of mars and saturn.
enum HelloGPCSP { jupiter, mars, saturn, venus, root };

// Our tree is
// (jupiter:0.113,(mars:0.15,saturn:0.1)venus:0.22):0.;
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
      EvolveRootward{PLV::ancestor_root, PLV::ancestor_leaf, HelloGPCSP::venus},
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

TEST_CASE("GPInstance: two pass optimization") {
  // Here t is the rootsplit (jupiter, other), and s is the osubsplit (mars, jupiter).
  enum PLV {
    p_jupiter,     // p at fake subsplit with jupiter
    p_mars,        // p at fake subsplit with mars
    p_saturn,      // p at fake subsplit with mars
    phat_stilde,   // evolved message coming from mars
    phat_s,        // evolved message coming from saturn
    p_s,           // p at venus
    phat_ttilde,   // evolved message coming from jupiter
    phat_t,        // evolved message coming from venus
    p_t,           // p at root
    rhat_t,        // stationary distribution "beyond" root
    r_ttilde,      // PLV pointing towards jupiter from root
    r_t,           // PLV pointing towards venus from root
    rhat_s,        // evolved message coming toward venus
    r_stilde,      // PLV pointing towards mars from venus
    r_s,           // PLV pointing towards saturn from venus
    rhat_jupiter,  // evolved message coming towards jupiter
    rhat_mars,     // evolved message coming towards mars
    rhat_saturn,   // evolved message coming towards saturn
  };

  GPOperationVector two_pass_likelihood_computation{
      // Evolve jupiter rootward.
      EvolveRootward{PLV::phat_ttilde, PLV::p_jupiter, HelloGPCSP::jupiter},
      // Evolve mars rootward.
      EvolveRootward{PLV::phat_stilde, PLV::p_mars, HelloGPCSP::mars},
      // Evolve saturn rootward.
      EvolveRootward{PLV::phat_s, PLV::p_saturn, HelloGPCSP::saturn},
      // Get rootward PLV for venus.
      Multiply{PLV::p_s, PLV::phat_stilde, PLV::phat_s},
      // Evolve venus rootward.
      EvolveRootward{PLV::phat_t, PLV::p_s, HelloGPCSP::venus},
      // Get rootward PLV for root.
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
      // Get message going towards venus.
      Multiply{PLV::r_t, PLV::phat_ttilde, PLV::rhat_t},
      // Evolve the message towards venus.
      EvolveLeafward{PLV::rhat_s, PLV::r_t, HelloGPCSP::venus},
      // Calculate likelihood at venus.
      Likelihood{HelloGPCSP::venus, PLV::rhat_s, PLV::p_s},
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

  GPOperationVector two_pass_optimization{
      // Optimize jupiter branch length and assign PLV on the root side.
      OptimizeRootward{PLV::phat_ttilde, PLV::p_jupiter, PLV::r_ttilde,
                       HelloGPCSP::jupiter},
      // Optimize mars branch length and assign PLV on the root side.
      OptimizeRootward{PLV::phat_stilde, PLV::p_mars, PLV::r_stilde, HelloGPCSP::mars},
      // Optimize saturn branch length and assign PLV on the root side.
      OptimizeRootward{PLV::phat_s, PLV::p_saturn, PLV::r_s, HelloGPCSP::saturn},
      // Get rootward PLV for venus.
      Multiply{PLV::p_s, PLV::phat_stilde, PLV::phat_s},
      // Optimize venus branch length and assign PLV on the root side.
      OptimizeRootward{PLV::phat_t, PLV::p_s, PLV::r_t, HelloGPCSP::venus},
      // Get rootward PLV for root.
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
      // Get message going towards venus.
      Multiply{PLV::r_t, PLV::phat_ttilde, PLV::rhat_t},
      // Evolve the message towards venus.
      EvolveLeafward{PLV::rhat_s, PLV::r_t, HelloGPCSP::venus},
      // Calculate likelihood at venus.
      Likelihood{HelloGPCSP::venus, PLV::rhat_s, PLV::p_s},
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

  // Test that the two-pass likelihood computation does the right thing at each internal
  // node.
  engine->ProcessOperations(two_pass_likelihood_computation);
  for (size_t idx = 0; idx <= root; idx++) {
    CHECK_LT(fabs(engine->GetLogLikelihoods()(idx) - -84.77961943), 1e-6);
  }

  // Test of our log likelihood derivative code on the jupiter branch.
  auto jupiter_optimization = OptimizeRootward{PLV::phat_ttilde, PLV::p_jupiter,
                                               PLV::r_ttilde, HelloGPCSP::jupiter};
  EigenVectorXd branch_lengths = engine->GetBranchLengths();
  auto [original_log_likelihood, log_likelihood_derivative] =
      engine->LogLikelihoodAndDerivative(jupiter_optimization);
  double branch_length_difference = 1e-7;
  branch_lengths(0) += branch_length_difference;
  engine->SetBranchLengths(branch_lengths);
  engine->ProcessOperations(two_pass_likelihood_computation);
  auto log_likelihood = engine->GetLogLikelihoods()(0);
  auto derivative_estimate =
      (log_likelihood - original_log_likelihood) / branch_length_difference;
  CHECK_LT(fabs(derivative_estimate - log_likelihood_derivative), 1e-6);

  // Trying out optimization.
  for (size_t pass_idx = 0; pass_idx < 8; ++pass_idx) {
    engine->ProcessOperations(two_pass_optimization);
  }
  std::cout << engine->GetBranchLengths() << std::endl;
  std::cout << engine->GetLogLikelihoods() << std::endl;
}
