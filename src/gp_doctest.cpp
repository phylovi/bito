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
  branch_lengths << 0, 0.22, 0.113, 0.15, 0.1;
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
  branch_lengths << 0, 0.22, 0.113, 0.15, 0.1;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  return inst;
}

GPInstance MakeHelloGPInstanceTwoTrees() {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/hello.fasta");
  inst.ReadNewickFile("data/hello_rooted_two_trees.nwk");
  inst.MakeEngine();
  inst.GetEngine()->SetBranchLengthsToConstant(1.);
  return inst;
}

EigenVectorXd MakeHelloGPInstanceOptimalBranchLengths() {
  EigenVectorXd hello_gp_optimal_branch_lengths(10);
  hello_gp_optimal_branch_lengths << 1, 1, 0.0540515854, 0.0150759956, 0.0158180779,
      0.0736678167, 0.206134746, 0.206803557, 0.0736678167, 0.0540515854;

  return hello_gp_optimal_branch_lengths;
}

// Hotstart
// Might need something like this..
GPInstance MakeHotStartGPInstance(){
    GPInstance inst("_ignore/mmapped_plv.data");
    inst.ReadFastaFile("data/hotstart.fasta");
    inst.ReadNewickFile("data/hotstart_bootstrap_sample.nwk");
    inst.MakeEngine();
    inst.PrintStatus();
    EigenVectorXd branch_lengths = inst.GetEngine()->GetBranchLengths();
    std::cout << "Here is the vector " << branch_lengths << std::endl;
    return inst;
}

// Outputting the test case branch length vector
EigenVectorXd HotStartTestCase(){
    Driver driver;
    RootedTreeCollection tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile("data/hotstart_bootstrap_sample.nwk"));

    const auto tree_count = tree_collection_.TreeCount();
    EigenVectorXd hotstart_test_case(tree_count);
    size_t i = 0;

    for (const auto& tree : tree_collection_.Trees()){
        hotstart_test_case(i) = tree.branch_lengths_[0];
        i++;
    }
    return hotstart_test_case;
}

TEST_CASE("GPInstance: straightforward classical likelihood calculation") {
  auto inst = MakeHelloGPInstance();
  auto engine = inst.GetEngine();

  inst.PopulatePLVs();
  inst.PrintDAG();
  inst.ComputeLikelihoods();

  EigenVectorXd realized_log_likelihoods = inst.GetEngine()->GetLogLikelihoods();
  CheckVectorXdEquality(-84.77961943, realized_log_likelihoods, 1e-6);

  CHECK_LT(fabs(engine->GetLogMarginalLikelihood() - -84.77961943), 1e-6);
}

TEST_CASE("GPInstance: marginal likelihood calculation") {
  auto inst = MakeHelloGPInstanceTwoTrees();
  auto engine = inst.GetEngine();

  auto branch_lengths = MakeHelloGPInstanceOptimalBranchLengths();
  engine->SetBranchLengths(branch_lengths);

  inst.PopulatePLVs();
  inst.ComputeLikelihoods();

  CHECK_LT(fabs(engine->GetLogMarginalLikelihood() - -80.6906345), 1e-6);
}

TEST_CASE("GPInstance: gradient calculation") {
  auto inst = MakeHelloGPInstanceSingleNucleotide();
  auto engine = inst.GetEngine();

  inst.PopulatePLVs();
  inst.ComputeLikelihoods();

  size_t root_idx = root;
  size_t child_idx = jupiter;
  size_t hello_node_count = 5;
  size_t root_jupiter_idx = 2;

  size_t leafward_idx =
      GPDAG::GetPLVIndexStatic(GPDAG::PLVType::P, hello_node_count, child_idx);
  size_t rootward_idx =
      GPDAG::GetPLVIndexStatic(GPDAG::PLVType::R, hello_node_count, root_idx);
  OptimizeBranchLength op{leafward_idx, rootward_idx, root_jupiter_idx};
  DoublePair log_lik_and_derivative = engine->LogLikelihoodAndDerivative(op);
  // Expect log lik: -4.806671945.
  // Expect log lik derivative: -0.6109379521.
  CHECK_LT(fabs(log_lik_and_derivative.first - -4.806671945), 1e-6);
  CHECK_LT(fabs(log_lik_and_derivative.second - -0.6109379521), 1e-6);
}

TEST_CASE("GPInstance: branch length optimization") {
  auto inst = MakeHelloGPInstanceTwoTrees();

  EigenVectorXd expected_branch_lengths = MakeHelloGPInstanceOptimalBranchLengths();
  inst.GetEngine()->SetBranchLengths(expected_branch_lengths);
  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  double expected_log_marginal = inst.GetEngine()->GetLogMarginalLikelihood();

  // Reset.
  inst = MakeHelloGPInstanceTwoTrees();
  inst.EstimateBranchLengths(1e-6, 100);
  EigenVectorXd realized_branch_lengths = inst.GetEngine()->GetBranchLengths();
  CheckVectorXdEquality(expected_branch_lengths, realized_branch_lengths, 1e-6);
  CHECK_LT(fabs(expected_log_marginal - -80.6906345), 1e-6);
}

GPInstance MakeFluAGPInstance(double rescaling_threshold) {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/fluA.fa");
  inst.ReadNewickFile("data/fluA.tree");
  inst.MakeEngine(rescaling_threshold);
  inst.GetEngine()->SetBranchLengthsToConstant(0.01);
  return inst;
}

double MakeAndRunFluAGPInstance(double rescaling_threshold) {
  auto inst = MakeFluAGPInstance(rescaling_threshold);
  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  return inst.GetEngine()->GetLogMarginalLikelihood();
}

TEST_CASE("GPInstance: rescaling") {
  double difference = MakeAndRunFluAGPInstance(GPEngine::default_rescaling_threshold_) -
                      MakeAndRunFluAGPInstance(1e-4);
  CHECK_LT(fabs(difference), 1e-10);
}

GPInstance MakeFiveTaxonRootedInstance() {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/five_taxon_rooted.fasta");
  inst.ReadNewickFile("data/five_taxon_rooted.nwk");
  inst.MakeEngine();
  return inst;
}

TEST_CASE("GPInstance: generate all trees") {
  auto inst = MakeFiveTaxonRootedInstance();
  auto rooted_tree_collection = inst.GenerateCompleteRootedTreeCollection();
  CHECK_EQ(rooted_tree_collection.TreeCount(), 4);
  CHECK_EQ(rooted_tree_collection.TopologyCounter().size(), 4);
}

EigenVectorXd MakeHotStartExpectedBranchLengths() {
  EigenVectorXd hotstart_expected_branch_lengths(8);
  hotstart_expected_branch_lengths << 0.135747939, 0.091601424, 0.11377897, 0.081359048,
    0.115477758, 0.081239333, 0.126825273, 0.135747939;

  return hotstart_expected_branch_lengths;
}

TEST_CASE("GPInstance: Hotstart branch lengths"){
    EigenVectorXd hotstart_test_case = HotStartTestCase();
    double true_mean = hotstart_test_case.array().mean();

    // TODO Figure out why hotstart function is only outputting 1's.
    // then it should be straightforward to finish the unit test
}
