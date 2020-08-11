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

TEST_CASE("GPInstance: Hotstart branch lengths") {
  const std::string tree_path = "data/hotstart_bootstrap_sample.nwk";
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/hotstart.fasta");
  inst.ReadNewickFile(tree_path);
  inst.MakeEngine();

  Driver driver;
  RootedTreeCollection tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(tree_path));

  const auto tree_count = tree_collection_.TreeCount();

  // » nw_topology data/hotstart_bootstrap_sample.nwk | nw_order - | sort | uniq -c
  // 1 (outgroup,(((z0,z1),z2),z3));
  // 33 (outgroup,((z0,z1),(z2,z3)));

  // This below is the GPCSP with sister z2, z3 and children z0, z1
  // That corresponds to indices 1, 2, and children 4, 3
  // 0110000011|0001000001, 2
  // [(1, z2), (2, z3), (3, z1), (4, z0), (0, outgroup)]
  // TODO a test that shows that in this case the IdTaxonMap is laid out like this.
  // Compare tree_collection_.IdTaxonMap is exactly this.

  EigenVectorXd hotstart_expected_branch_lengths(33);
  hotstart_expected_branch_lengths << 0.1175370000, 0.1175750000, 0.1195780000,
      0.0918962000, 0.0918931000, 0.1192590000, 0.0906988000, 0.0906972000,
      0.0905154000, 0.0903663000, 0.1245620000, 0.1244890000, 0.1245050000,
      0.1245550000, 0.1245680000, 0.1248920000, 0.1248490000, 0.1164070000,
      0.1164110000, 0.1164120000, 0.1245670000, 0.1245650000, 0.1245670000,
      0.1245670000, 0.1240790000, 0.1242540000, 0.1242160000, 0.1242560000,
      0.1892030000, 0.1894900000, 0.1895430000, 0.1896900000, 0.1905710000;

  double true_mean = hotstart_expected_branch_lengths.array().mean();
  std::cout << "true mean " << true_mean << std::endl;
  inst.PrintGPCSPIndexer();

  inst.HotStartBranchLengths();
  std::cout << "After hot start: " << inst.GetEngine()->GetBranchLengths() << std::endl;

  CHECK_EQ(true_mean, inst.GetEngine()->GetBranchLengths()(2));
}
