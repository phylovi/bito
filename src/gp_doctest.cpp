#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "gp_instance.hpp"

using namespace GPOperations;

// GPCSP stands for generalized PCSP-- see text.
// Let the "venus" node be the common ancestor of mars and saturn.
enum HelloGPCSP { jupiter, mars, saturn, venus, root };

// Our tree is
// (jupiter:0.113,(mars:0.15,saturn:0.1)venus:0.22):0.;
// You can see a helpful diagram at
// https://github.com/phylovi/libsbn/issues/213#issuecomment-624195267
GPInstance MakeHelloGPInstance() {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/hello.fasta");
  inst.ReadNewickFile("data/hello_rooted.nwk");
  inst.MakeEngine();
  EigenVectorXd branch_lengths(5);
  // Order set by HelloGPCSP.
  branch_lengths << 0, 0.1, 0.15, 0.113, 0.22;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  return inst;
}

GPInstance MakeHelloGPInstanceSingleNucleotide() {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/hello_single_nucleotide.fasta");
  inst.ReadNewickFile("data/hello_rooted.nwk");
  inst.MakeEngine();
  EigenVectorXd branch_lengths(5);
  // Order set by HelloGPCSP.
  branch_lengths << 0, 0.1, 0.15, 0.113, 0.22;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  return inst;
}

GPInstance MakeHelloGPInstanceTwoTrees() {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/hello.fasta");
  inst.ReadNewickFile("data/hello_rooted_two_trees.nwk");
  inst.MakeEngine();
  EigenVectorXd branch_lengths =
      EigenVectorXd::Ones(inst.GetEngine()->GetBranchLengths().size());
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  return inst;
}

TEST_CASE("GPInstance: straightforward classical likelihood calculation") {
  auto inst = MakeHelloGPInstance();
  // auto engine = inst.GetEngine();

  // inst.PopulatePLVs();
  // inst.PrintDAG();
  // inst.ComputeLikelihoods();

  // CHECK_LT(fabs(engine->GetLogMarginalLikelihood() - -84.77961943), 1e-6);
}

// TEST_CASE("GPInstance: marginal likelihood calculation") {
//   auto inst = MakeHelloGPInstanceTwoTrees();
//   auto engine = inst.GetEngine();
//
//   inst.PopulatePLVs();
//   inst.ComputeLikelihoods();
//   std::cout << engine->GetLogMarginalLikelihood() << std::endl;
//   inst.EstimateBranchLengths(1e-6, 100);
//
//   CHECK_LT(fabs(engine->GetLogMarginalLikelihood() - -79.9944001), 1e-6);
// }

// TEST_CASE("GPInstance: gradient calculation") {
//   auto inst = MakeHelloGPInstanceSingleNucleotide();
//   auto engine = inst.GetEngine();
//
//   inst.PopulatePLVs();
//   inst.ComputeLikelihoods();
//
//   size_t root_idx = root;
//   size_t child_idx = jupiter;
//   size_t hello_node_count = 5;
//   size_t leafward_idx = GetPLVIndex(PLVType::P, hello_node_count, child_idx);
//   size_t rootward_idx = GetPLVIndex(PLVType::R, hello_node_count, root_idx);
//   size_t pcsp_idx = inst.GetPCSPIndex(root_idx, child_idx, false);
//   OptimizeBranchLength op{leafward_idx, rootward_idx, pcsp_idx};
//   DoublePair log_lik_and_derivative = engine->LogLikelihoodAndDerivative(op);
//   // Expect log lik: -4.806671945.
//   // Expect log lik derivative: -0.6109379521.
//   CHECK_LT(fabs(log_lik_and_derivative.first - -4.806671945), 1e-6);
//   CHECK_LT(fabs(log_lik_and_derivative.second - -0.6109379521), 1e-6);
// }
