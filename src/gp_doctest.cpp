// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

// ** Doctest include must go first for all header tests to run.
#include "doctest.h"
// **

#include "combinatorics.hpp"
#include "gp_instance.hpp"
#include "phylo_model.hpp"
#include "unrooted_sbn_instance.hpp"

using namespace GPOperations;  // NOLINT

// GPCSP stands for generalized PCSP-- see text.

// Let the "venus" node be the common ancestor of mars and saturn.
enum HelloGPCSP { jupiter, mars, saturn, venus, root };

// Compute the exact marginal likelihood via brute force to compare with generalized
// pruning. We assume that the trees in `newick_path` are all of the trees over which we
// should marginalize.
double ComputeExactMarginal(const std::string& newick_path,
                            const std::string& fasta_path,
                            bool with_tree_prior = true) {
  UnrootedSBNInstance sbn_instance("charlie");
  sbn_instance.ReadNewickFile(newick_path);
  size_t tree_count = sbn_instance.TreeCount();
  const Alignment alignment = Alignment::ReadFasta(fasta_path);
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  sbn_instance.SetAlignment(alignment);
  sbn_instance.PrepareForPhyloLikelihood(simple_specification, 1);

  double exact_marginal_log_lik = 0.0;
  for (size_t i = 0; i < alignment.Length(); i++) {
    sbn_instance.SetAlignment(alignment.ExtractSingleColumnAlignment(i));
    sbn_instance.PrepareForPhyloLikelihood(simple_specification, 1);
    double log_val = DOUBLE_NEG_INF;
    for (double val : sbn_instance.LogLikelihoods()) {
      log_val = NumericalUtils::LogAdd(log_val, val);
    }
    log_val += with_tree_prior ? log(1. / tree_count) : 0.0;
    exact_marginal_log_lik += log_val;
  }
  return exact_marginal_log_lik;
}

GPInstance GPInstanceOfFiles(const std::string& fasta_path,
                             const std::string& newick_path) {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile(fasta_path);
  inst.ReadNewickFile(newick_path);
  inst.MakeEngine();
  return inst;
}

// Our tree is (see check below)
// (jupiter:0.113,(mars:0.15,saturn:0.1)venus:0.22):0.;
// You can see a helpful diagram at
// https://github.com/phylovi/libsbn/issues/213#issuecomment-624195267
GPInstance MakeHelloGPInstance(const std::string& fasta_path) {
  auto inst = GPInstanceOfFiles(fasta_path, "data/hello_rooted.nwk");
  EigenVectorXd branch_lengths(5);
  // Order set by HelloGPCSP.
  branch_lengths << 0, 0.22, 0.113, 0.15, 0.1;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  CHECK_EQ(inst.GenerateCompleteRootedTreeCollection().Newick(),
           "(jupiter:0.113,(mars:0.15,saturn:0.1):0.22):0;\n");
  return inst;
}

GPInstance MakeHelloGPInstance() { return MakeHelloGPInstance("data/hello.fasta"); }

GPInstance MakeHelloGPInstanceSingleNucleotide() {
  return MakeHelloGPInstance("data/hello_single_nucleotide.fasta");
}

GPInstance MakeHelloGPInstanceTwoTrees() {
  auto inst = GPInstanceOfFiles("data/hello.fasta", "data/hello_rooted_two_trees.nwk");
  inst.GetEngine()->SetBranchLengthsToConstant(1.);
  return inst;
}

// The sequences for this were obtained by cutting DS1 down to 5 taxa by taking the
// first 4 taxa then moving taxon 15 (Latimera) to be number 5. The alignment was
// trimmed to 500 sites by using seqmagick convert with `--cut 500:1000`.
GPInstance MakeDS1Reduced5Instance() {
  auto inst = GPInstanceOfFiles("data/ds1-reduced-5.fasta", "data/ds1-reduced-5.nwk");
  return inst;
}

EigenVectorXd MakeHelloGPInstanceMarginalLikelihoodTestBranchLengths() {
  EigenVectorXd hello_gp_optimal_branch_lengths(10);
  hello_gp_optimal_branch_lengths << 1, 1, 0.066509261, 0.00119570257, 0.00326456973,
      0.0671995398, 0.203893516, 0.204056242, 0.0669969961, 0.068359082;

  return hello_gp_optimal_branch_lengths;
}

EigenVectorXd MakeHelloGPInstanceRegressionTestBranchLengths() {
  EigenVectorXd hello_gp_optimal_branch_lengths(10);
  hello_gp_optimal_branch_lengths << 1, 1, 0.066509261, 0.00119570257, 0.00326456973,
      0.0671995398, 0.203893516, 0.204056242, 0.0669969961, 0.068359082;

  return hello_gp_optimal_branch_lengths;
}

void CheckBranchLengthOptimization(GPInstance& inst) {
  auto old_operations = inst.GetDAG().BranchLengthOptimization();
  auto new_operations = inst.GetDAG().NewBranchLengthOptimization();
  CHECK_EQ(GenericToString(old_operations), GenericToString(new_operations));
}

TEST_CASE("GPInstance: straightforward classical likelihood calculation") {
  auto inst = MakeHelloGPInstance();
  auto engine = inst.GetEngine();

  inst.ResetMarginalLikelihoodAndPopulatePLVs();
  inst.ComputeLikelihoods();

  EigenVectorXd realized_log_likelihoods =
      inst.GetEngine()->GetPerGPCSPLogLikelihoods();
  CheckVectorXdEquality(-84.77961943, realized_log_likelihoods, 1e-6);

  CHECK_LT(fabs(engine->GetLogMarginalLikelihood() - -84.77961943), 1e-6);

  CheckBranchLengthOptimization(inst);
}

TEST_CASE("GPInstance: two tree marginal likelihood calculation") {
  auto inst = MakeHelloGPInstanceTwoTrees();
  auto engine = inst.GetEngine();

  auto branch_lengths = MakeHelloGPInstanceMarginalLikelihoodTestBranchLengths();
  engine->SetBranchLengths(branch_lengths);

  inst.ResetMarginalLikelihoodAndPopulatePLVs();
  inst.ComputeLikelihoods();
  double gp_marginal_log_likelihood = engine->GetLogMarginalLikelihood();
  double exact_log_likelihood =
      ComputeExactMarginal("data/hello_unrooted_two_trees.nwk", "data/hello.fasta");
  CHECK_LT(fabs(gp_marginal_log_likelihood - exact_log_likelihood), 1e-6);

  CheckBranchLengthOptimization(inst);
}

TEST_CASE("GPInstance: DS1-reduced-5 marginal likelihood calculation") {
  auto inst = MakeDS1Reduced5Instance();
  auto engine = inst.GetEngine();

  inst.EstimateBranchLengths(0.0001, 100, true);
  inst.ResetMarginalLikelihoodAndPopulatePLVs();
  inst.ComputeLikelihoods();
  double gp_marginal_log_likelihood = engine->GetLogMarginalLikelihood();

  const auto trees = inst.CurrentlyLoadedTreesWithGPBranchLengths();
  // "ds1-reduced-5.gp-branch-lengths.unrooted.nwk" created by
  // trees.ToNewickFile("data/ds1-reduced-5.gp-branch-lengths.nwk");
  // then derooting with `nw_reroot -d`.
  double exact_log_likelihood = ComputeExactMarginal(
      "data/ds1-reduced-5.gp-branch-lengths.unrooted.nwk", "data/ds1-reduced-5.fasta");
  CHECK_LT(fabs(gp_marginal_log_likelihood - exact_log_likelihood), 1e-6);

  CheckBranchLengthOptimization(inst);
}

TEST_CASE("GPInstance: gradient calculation") {
  auto inst = MakeHelloGPInstanceSingleNucleotide();
  auto engine = inst.GetEngine();

  inst.ResetMarginalLikelihoodAndPopulatePLVs();
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
  inst.ResetMarginalLikelihoodAndPopulatePLVs();
  inst.ComputeLikelihoods();
  return inst.GetEngine()->GetLogMarginalLikelihood();
}

// Regression test.
TEST_CASE("GPInstance: branch length optimization") {
  auto inst = MakeHelloGPInstanceTwoTrees();

  EigenVectorXd expected_branch_lengths =
      MakeHelloGPInstanceRegressionTestBranchLengths();

  // Reset.
  inst = MakeHelloGPInstanceTwoTrees();
  inst.EstimateBranchLengths(1e-6, 100, true);
  EigenVectorXd realized_branch_lengths = inst.GetEngine()->GetBranchLengths();
  CheckVectorXdEquality(expected_branch_lengths, realized_branch_lengths, 1e-6);
}

TEST_CASE("GPInstance: rescaling") {
  double difference = MakeAndRunFluAGPInstance(GPEngine::default_rescaling_threshold_) -
                      MakeAndRunFluAGPInstance(1e-4);
  CHECK_LT(fabs(difference), 1e-10);
}

TEST_CASE("GPInstance: hotstart branch lengths") {
  // » nw_topology data/hotstart_bootstrap_sample.nwk | nw_order - | sort | uniq -c
  // 1 (outgroup,(((z0,z1),z2),z3));
  // 33 (outgroup,((z0,z1),(z2,z3)));
  const std::string tree_path = "data/hotstart_bootstrap_sample.nwk";
  GPInstance inst("_ignore/mmapped_plv.data");
  // This is just a dummy fasta file, which is required to make an Engine.
  inst.ReadFastaFile("data/hotstart.fasta");
  inst.ReadNewickFile(tree_path);
  inst.MakeEngine();

  // We are going to verify correct assignment of the PCSP with sister z2, z3 and
  // children z0, z1, which only appears in the tree (outgroup,((z0,z1),(z2,z3))).
  // Vector of taxon names: [outgroup, z2, z3, z1, z0]
  // So, this below is the desired GPCSP (in full subsplit notation), which corresponds
  // to sister indices 1, 2, and children 4, 3:
  // 0110000011|0001000001, 2
  // Thus we are interested in the branch length index 2.

  // These branch lengths are obtained by excluding (outgroup,(((z0,z1),z2),z3)) (which
  // doesn't have this PCSP) and grabbing the rest of the branch lengths.
  EigenVectorXd hotstart_expected_branch_lengths(33);
  hotstart_expected_branch_lengths << 0.1175370000, 0.1175750000, 0.1195780000,
      0.0918962000, 0.0918931000, 0.1192590000, 0.0906988000, 0.0906972000,
      0.0905154000, 0.0903663000, 0.1245620000, 0.1244890000, 0.1245050000,
      0.1245550000, 0.1245680000, 0.1248920000, 0.1248490000, 0.1164070000,
      0.1164110000, 0.1164120000, 0.1245670000, 0.1245650000, 0.1245670000,
      0.1245670000, 0.1240790000, 0.1242540000, 0.1242160000, 0.1242560000,
      0.1892030000, 0.1894900000, 0.1895430000, 0.1896900000, 0.1905710000;
  double true_mean = hotstart_expected_branch_lengths.array().mean();
  inst.HotStartBranchLengths();
  CHECK_EQ(true_mean, inst.GetEngine()->GetBranchLengths()(2));
}

GPInstance MakeFiveTaxaInstance() {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/five_taxon.fasta");
  inst.ReadNewickFile("data/five_taxon_rooted.nwk");
  inst.MakeEngine();
  return inst;
}

TEST_CASE("GPInstance: generate all trees") {
  auto inst = MakeFiveTaxaInstance();
  auto rooted_tree_collection = inst.GenerateCompleteRootedTreeCollection();
  CHECK_EQ(rooted_tree_collection.TreeCount(), 4);
  CHECK_EQ(rooted_tree_collection.TopologyCounter().size(), 4);
}

TEST_CASE("GPInstance: marginal likelihood on five taxa") {
  auto inst = MakeFiveTaxaInstance();
  inst.EstimateBranchLengths(1e-6, 10, true);
  double gp_marginal_log_likelihood = inst.GetEngine()->GetLogMarginalLikelihood();
  double exact_marginal_log_likelihood = ComputeExactMarginal(
      "data/five_taxon_unrooted_with_branch_lengths.nwk", "data/five_taxon.fasta");
  CHECK_LT(abs(exact_marginal_log_likelihood - gp_marginal_log_likelihood), 1e-6);
}

TEST_CASE("GPInstance: test populate PLV") {
  // This test makes sure that ResetMarginalLikelihoodAndPopulatePLVs correctly
  // re-populates the PLVs using the current branch lengths.
  auto inst = MakeFiveTaxaInstance();
  inst.EstimateBranchLengths(1e-6, 10, true);
  inst.ComputeLikelihoods();
  size_t length = inst.GetEngine()->GetLogLikelihoodMatrix().rows();
  const EigenVectorXd log_likelihoods1 =
      inst.GetEngine()->GetPerGPCSPLogLikelihoods(0, length);
  inst.ResetMarginalLikelihoodAndPopulatePLVs();
  inst.ComputeLikelihoods();
  const EigenVectorXd log_likelihoods2 = inst.GetEngine()->GetPerGPCSPLogLikelihoods();
  CheckVectorXdEquality(log_likelihoods1, log_likelihoods2, 1e-6);
}

TEST_CASE("GPInstance: SBN root split probabilities on five taxa") {
  auto inst = MakeFiveTaxaInstance();
  inst.GetEngine()->SetBranchLengthsToConstant(0.1);
  inst.ResetMarginalLikelihoodAndPopulatePLVs();
  // We need to call ComputeLikelihoods to populate the likelihood matrix.
  // Note: EstimateBranchLengths doesn't populate the likelihood matrix.
  inst.ComputeLikelihoods();

  EigenVectorXd log_likelihood_vector = inst.GetEngine()->GetPerGPCSPLogLikelihoods();

  // Let s be a subsplit and k be the site. Then,
  // log_likelihood_matrix.row(s)[k] =
  //    \log \sum_{\tau : s \in \tau} q(\tau) P(y_k | \tau),
  // log_likelihood_vector[s] =
  //    \sum_{k=1}^{K} \log \sum_{\tau : s \in \tau} q(\tau) P(y_k | \tau).
  // To test this, we are going to compute P(y_k | \tau) for {\tau : s \in \tau} and
  // multiply this by q(\tau) = 1/4 since we are assuming uniform prior.

  // The collection of trees that we are looking at has 3 rootplits where one root split
  // generates two trees and the other 2 root splits generating one tree each
  // for the total of 4 trees.

  // We will compare the values against the 3 rootsplits, since we cannot assume
  // the ordering due to different implementation of the map, we will sort the values
  // before comparison.

  double log_lik_tree_1 =
      ComputeExactMarginal("data/five_taxon_tree1.nwk", "data/five_taxon.fasta", false);
  double log_lik_tree_2 =
      ComputeExactMarginal("data/five_taxon_tree2.nwk", "data/five_taxon.fasta", false);
  double log_lik_trees_3_4 = ComputeExactMarginal("data/five_taxon_trees_3_4.nwk",
                                                  "data/five_taxon.fasta", false);

  size_t alignment_length = 4;
  // Uniform prior over 4 trees.
  double log_q = log(1. / 4);
  EigenVectorXd expected_log_lik_vector_at_rootsplits(3);
  expected_log_lik_vector_at_rootsplits[0] = log_lik_tree_1 + alignment_length * log_q;
  expected_log_lik_vector_at_rootsplits[1] = log_lik_tree_2 + alignment_length * log_q;
  expected_log_lik_vector_at_rootsplits[2] =
      log_lik_trees_3_4 + alignment_length * log_q;
  EigenVectorXd realized_log_lik_vector_at_rootsplits =
      log_likelihood_vector.segment(0, 3);
  CheckVectorXdEqualityAfterSorting(realized_log_lik_vector_at_rootsplits,
                                    expected_log_lik_vector_at_rootsplits, 1e-6);

  inst.EstimateSBNParameters();
  EigenVectorXd realized_q = inst.GetEngine()->GetSBNParameters().segment(0, 3);
  // The expected values for the SBN parameters: q[s] \propto log_lik[s] + log_prior[s].
  // The SBN params are initialized so that we get a uniform distribution over the
  // trees. For the rootsplits, the values are (1/4, 1/4, 2/4) corresponding to the
  // entries in expected_log_lik_vector_at_rootsplits.
  EigenVectorXd log_prior(3);
  log_prior << log(1. / 4), log(1. / 4), log(2. / 4);
  EigenVectorXd expected_q = expected_log_lik_vector_at_rootsplits + log_prior;
  NumericalUtils::ProbabilityNormalizeInLog(expected_q);
  expected_q = expected_q.array().exp();
  CheckVectorXdEqualityAfterSorting(realized_q, expected_q, 1e-6);
}

TEST_CASE("GPInstance: CurrentlyLoadedTreesWithGPBranchLengths") {
  auto inst = MakeHelloGPInstanceSingleNucleotide();
  EigenVectorXd branch_lengths(5);
  branch_lengths << 0, 0.1, 0.2, 0.3, 0.4;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  auto trees = inst.CurrentlyLoadedTreesWithGPBranchLengths();
  CHECK_EQ(trees.Newick(), "(jupiter:0.2,(mars:0.3,saturn:0.4):0.1):0;\n");
}

TEST_CASE("GPInstance: CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths") {
  GPInstance inst("_ignore/mmapped_plv.data");
  inst.ReadFastaFile("data/five_taxon.fasta");
  inst.ReadNewickFile("data/five_taxon_rooted_more.nwk");
  inst.MakeEngine();
  inst.GetEngine()->SetBranchLengthsToConstant(0.9);
  // Only take trees that have (x4,(x2,x3)).
  auto trees =
      inst.CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths("000010011000010");
  CHECK_EQ(trees.Newick(),
           "((x0:0.9,x1:0.9):0.9,((x2:0.9,x3:0.9):0.9,x4:0.9):0.9):0;\n"
           "(x0:0.9,(x1:0.9,((x2:0.9,x3:0.9):0.9,x4:0.9):0.9):0.9):0;\n");
}

TEST_CASE("GPInstance: Priors") {
  auto inst = GPInstanceOfFiles("data/four-numbered-taxa.fasta",
                                "data/four-taxon-two-tree-rootsplit-uncertainty.nwk");
  // Here are the trees:
  // (((1,2),3),4);
  // ((1,(2,3)),4);
  // ((1,2),(3,4));
  //
  // Here's the interesting part of the indexer:
  // 0001|1110,      0
  // 0011|1100,      1
  // 0001|1110|0110, 4
  // 0001|1110|0010, 5
  auto support = inst.GetDAG().BuildUniformOnTopologicalSupportPrior();
  CHECK_LT(fabs(support[0] - 2. / 3.), 1e-10);
  CHECK_LT(fabs(support[1] - 1. / 3.), 1e-10);
  CHECK_LT(fabs(support[4] - 1. / 2.), 1e-10);
  CHECK_LT(fabs(support[5] - 1. / 2.), 1e-10);
  auto all = inst.GetDAG().BuildUniformOnAllTopologiesPrior();
  // There are 15 topologies on 4 taxa.
  // There are 3 topologies on 3 taxa, so there are 3 topologies with rootsplit
  // 0001|1110.
  CHECK_LT(fabs(all[0] - 3. / 15.), 1e-10);
  // There is only 1 topology with rootsplit 0011|1100.
  CHECK_LT(fabs(all[1] - 1. / 15.), 1e-10);
  // There are 3 topologies on 3 taxa.
  CHECK_LT(fabs(all[4] - 1. / 3.), 1e-10);
  CHECK_LT(fabs(all[5] - 1. / 3.), 1e-10);
}
