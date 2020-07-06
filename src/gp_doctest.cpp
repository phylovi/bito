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
  branch_lengths << 0, 0.113, 0.22, 0.1, 0.15;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  return inst;
}

GPInstance MakeHelloGPInstance2() {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/hello.fasta");
  inst.ReadNewickFile("data/hello_rooted_two_trees.nwk");
  inst.MakeEngine();
  EigenVectorXd branch_lengths = EigenVectorXd::Ones(inst.GetEngine()->GetBranchLengths().size());
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  std::cout << inst.GetEngine()->GetBranchLengths() << std::endl;
  std::cout << inst.GetEngine()->GetSBNParameters() << std::endl;
  return inst;
}

TEST_CASE("GPInstance: straightforward classical likelihood calculation") {
  auto inst = MakeHelloGPInstance();
  auto engine = inst.GetEngine();

  inst.PopulatePLVs();
  inst.ComputeLikelihoods();

  CHECK_LT(fabs(engine->GetLogMarginalLikelihood() - -84.77961943), 1e-6);
}

TEST_CASE("GPInstance: marginal likelihood calculation") {
  auto inst = MakeHelloGPInstance2();
  auto engine = inst.GetEngine();

  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  std::cout << engine->GetLogMarginalLikelihood() << std::endl;
  inst.EstimateBranchLengths(1e-6, 10);

  CHECK_LT(fabs(engine->GetLogMarginalLikelihood() - -79.9944001), 1e-6);
}
