// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

// ** Doctest include must go first for all header tests to run.
#include "doctest.h"
// **

#include "combinatorics.hpp"
#include "gp_instance.hpp"
#include "phylo_model.hpp"
#include "rooted_sbn_instance.hpp"
#include "tidy_subsplit_dag.hpp"

using namespace GPOperations;  // NOLINT

// GPCSP stands for generalized PCSP-- see text.

// Let the "venus" node be the common ancestor of mars and saturn.
enum HelloGPCSP { jupiter, mars, saturn, venus, root };

// *** GPInstances used for testing ***

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
  return GPInstanceOfFiles("data/hello.fasta", "data/hello_rooted_two_trees.nwk");
}

GPInstance MakeFiveTaxonInstance() {
  return GPInstanceOfFiles("data/five_taxon.fasta", "data/five_taxon_rooted.nwk");
}

// The sequences for this were obtained by cutting DS1 down to 5 taxa by taking the
// first 4 taxa then moving taxon 15 (Latimera) to be number 5. The alignment was
// trimmed to 500 sites by using seqmagick convert with `--cut 500:1000`.
// The DAG obtained by `inst.SubsplitDAGToDot("_ignore/ds1-reduced-5.dot");` can be seen
// at
// https://user-images.githubusercontent.com/112708/106355238-8a0d7700-62ab-11eb-8830-f5e3cb6790e1.png
GPInstance MakeDS1Reduced5Instance() {
  auto inst = GPInstanceOfFiles("data/ds1-reduced-5.fasta", "data/ds1-reduced-5.nwk");
  return inst;
}

GPInstance MakeFluAGPInstance(double rescaling_threshold) {
  auto inst = GPInstanceOfFiles("data/fluA.fa", "data/fluA.tree");
  inst.MakeEngine(rescaling_threshold);
  inst.GetEngine()->SetBranchLengthsToConstant(0.01);
  return inst;
}

EigenVectorXd MakeHelloGPInstanceMarginalLikelihoodTestBranchLengths() {
  EigenVectorXd hello_gp_optimal_branch_lengths(10);
  hello_gp_optimal_branch_lengths << 1, 1, 0.066509261, 0.00119570257, 0.00326456973,
      0.0671995398, 0.203893516, 0.204056242, 0.0669969961, 0.068359082;

  return hello_gp_optimal_branch_lengths;
}

TEST_CASE("GPInstance: straightforward classical likelihood calculation") {
  auto inst = MakeHelloGPInstance();
  auto engine = inst.GetEngine();

  inst.PopulatePLVs();
  inst.ComputeLikelihoods();

  EigenVectorXd realized_log_likelihoods =
      inst.GetEngine()->GetPerGPCSPLogLikelihoods();
  CheckVectorXdEquality(-84.77961943, realized_log_likelihoods, 1e-6);

  CHECK_LT(fabs(engine->GetLogMarginalLikelihood() - -84.77961943), 1e-6);
}

// Compute the exact marginal likelihood via brute force to compare with generalized
// pruning.
// IMPORTANT: We assume that the trees in `newick_path` are all of the trees over which
// we should marginalize. So if you have generated a subsplit DAG with a set of trees,
// use GenerateCompleteRootedTreeCollection to get all the trees over which you will be
// marginalizing.
// If we rename things in #288, let's do that in the body of this function too.
std::pair<double, StringDoubleMap> ComputeExactMarginal(const std::string& newick_path,
                                                        const std::string& fasta_path) {
  RootedSBNInstance sbn_instance("charlie");
  sbn_instance.ReadNewickFile(newick_path);
  sbn_instance.ProcessLoadedTrees();
  const Alignment alignment = Alignment::ReadFasta(fasta_path);
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  sbn_instance.SetAlignment(alignment);
  sbn_instance.PrepareForPhyloLikelihood(simple_specification, 1);

  const size_t tree_count = sbn_instance.TreeCount();
  const size_t gpcsp_count = sbn_instance.SBNSupport().GPCSPCount();
  auto indexer_representations = sbn_instance.MakeIndexerRepresentations();

  double exact_marginal_log_lik = 0.0;
  EigenVectorXd exact_per_pcsp_log_marginals(gpcsp_count);
  exact_per_pcsp_log_marginals.setZero();
  double log_prior_term = log(1. / tree_count);

  for (size_t column_idx = 0; column_idx < alignment.Length(); column_idx++) {
    sbn_instance.SetAlignment(alignment.ExtractSingleColumnAlignment(column_idx));
    sbn_instance.PrepareForPhyloLikelihood(simple_specification, 1);
    auto per_site_phylo_likelihoods = sbn_instance.UnrootedLogLikelihoods();

    double per_site_log_marginal = DOUBLE_NEG_INF;
    EigenVectorXd per_site_per_pcsp_log_marginals(gpcsp_count);
    per_site_per_pcsp_log_marginals.setConstant(DOUBLE_NEG_INF);

    for (size_t tree_idx = 0; tree_idx < tree_count; tree_idx++) {
      const auto per_site_phylo_likelihood = per_site_phylo_likelihoods[tree_idx];
      per_site_log_marginal =
          NumericalUtils::LogAdd(per_site_log_marginal, per_site_phylo_likelihood);
      for (const auto& gpcsp_idx : indexer_representations.at(tree_idx)) {
        per_site_per_pcsp_log_marginals[gpcsp_idx] = NumericalUtils::LogAdd(
            per_site_per_pcsp_log_marginals[gpcsp_idx], per_site_phylo_likelihood);
      }
    }
    per_site_log_marginal += log_prior_term;
    per_site_per_pcsp_log_marginals.array() += log_prior_term;

    exact_marginal_log_lik += per_site_log_marginal;
    exact_per_pcsp_log_marginals.array() += per_site_per_pcsp_log_marginals.array();
  }
  return {
      exact_marginal_log_lik,
      UnorderedMapOf(sbn_instance.PrettyIndexedVector(exact_per_pcsp_log_marginals))};
}

void CheckExactMapVsGPVector(const StringDoubleMap& exact_map,
                             const StringDoubleVector& gp_vector) {
  for (const auto& [gp_string, gp_value] : gp_vector) {
    if (exact_map.find(gp_string) == exact_map.end()) {
      Assert(!Bitset(gp_string.substr(gp_string.rfind('|') + 1)).Any(),
             "Missing a non-fake bitset in CheckExactMapVsGPVector.");
    } else {
      const double tolerance = 1e-5;
      const double error = fabs(exact_map.at(gp_string) - gp_value);
      if (error > tolerance) {
        std::cout << "check failed for " << gp_string << ":" << std::endl;
      }
      CHECK_LT(error, tolerance);
    }
  }
}

// Test the composite marginal to that generated by ComputeExactMarginal.
//
// IMPORTANT: See the note about appropriate tree file input to that function, as the
// same applies here.
void TestCompositeMarginal(GPInstance inst, const std::string& fasta_path) {
  inst.EstimateBranchLengths(0.0001, 100, true);
  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  inst.ComputeMarginalLikelihood();
  std::string tree_path = "_ignore/test_marginal_trees.nwk";
  const auto trees = inst.CurrentlyLoadedTreesWithGPBranchLengths();
  trees.ToNewickFile(tree_path);

  auto [exact_log_likelihood, exact_per_pcsp_log_marginal] =
      ComputeExactMarginal(tree_path, fasta_path);
  double gp_marginal_log_likelihood = inst.GetEngine()->GetLogMarginalLikelihood();
  auto gp_per_pcsp_log_marginal =
      inst.PrettyIndexedPerGPCSPComponentsOfFullLogMarginal();
  CHECK_LT(fabs(gp_marginal_log_likelihood - exact_log_likelihood), 1e-6);
  CheckExactMapVsGPVector(exact_per_pcsp_log_marginal, gp_per_pcsp_log_marginal);
}

TEST_CASE("GPInstance: two tree marginal likelihood calculation") {
  TestCompositeMarginal(MakeHelloGPInstanceTwoTrees(), "data/hello.fasta");
}

TEST_CASE("GPInstance: marginal likelihood on five taxa") {
  TestCompositeMarginal(MakeFiveTaxonInstance(), "data/five_taxon.fasta");
}

TEST_CASE("GPInstance: DS1-reduced-5 marginal likelihood calculation") {
  TestCompositeMarginal(MakeDS1Reduced5Instance(), "data/ds1-reduced-5.fasta");
}

TEST_CASE("GPInstance: marginal likelihood on seven taxa and four trees") {
  const std::string fasta_path = "data/7-taxon-slice-of-ds1.fasta";
  // See the DAG at
  // https://github.com/phylovi/libsbn/issues/323#issuecomment-783349464
  TestCompositeMarginal(
      GPInstanceOfFiles(fasta_path, "data/simplest-hybrid-marginal-all-trees.nwk"),
      fasta_path);
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
  // TODO: The correct index on mac is 1, not 2. Need to make this test more general to pass on all systems.
  std::cout << "The hotstart branch lengths are: " << inst.GetEngine()->GetBranchLengths();
  CHECK_EQ(true_mean, inst.GetEngine()->GetBranchLengths()(2));
}

TEST_CASE("GPInstance: generate all trees") {
  auto inst = MakeFiveTaxonInstance();
  auto rooted_tree_collection = inst.GenerateCompleteRootedTreeCollection();
  CHECK_EQ(rooted_tree_collection.TreeCount(), 4);
  CHECK_EQ(rooted_tree_collection.TopologyCounter().size(), 4);
}

TEST_CASE("GPInstance: test populate PLV") {
  // This test makes sure that PopulatePLVs correctly
  // re-populates the PLVs using the current branch lengths.
  auto inst = MakeFiveTaxonInstance();
  inst.EstimateBranchLengths(1e-6, 10, true);
  inst.ComputeLikelihoods();
  size_t length = inst.GetEngine()->GetLogLikelihoodMatrix().rows();
  const EigenVectorXd log_likelihoods1 =
      inst.GetEngine()->GetPerGPCSPLogLikelihoods(0, length);
  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  const EigenVectorXd log_likelihoods2 = inst.GetEngine()->GetPerGPCSPLogLikelihoods();
  CheckVectorXdEquality(log_likelihoods1, log_likelihoods2, 1e-6);
}

TEST_CASE("GPInstance: SBN root split probabilities on five taxa") {
  auto inst = MakeFiveTaxonInstance();
  inst.GetEngine()->SetBranchLengthsToConstant(0.1);
  inst.PopulatePLVs();
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
  // multiply this by q(\tau) = 1/4 since we are assuming a uniform prior.

  // The collection of trees that we are looking at has 3 rootplits where one root split
  // generates two trees and the other 2 root splits generating one tree each
  // for the total of 4 trees.

  // We will compare the values against the 3 rootsplits, since we cannot assume
  // the ordering due to different implementation of the map, we will sort the values
  // before comparison.

  auto [log_lik_tree_1, ignored_1] =
      ComputeExactMarginal("data/five_taxon_tree1.nwk", "data/five_taxon.fasta");
  auto [log_lik_tree_2, ignored_2] =
      ComputeExactMarginal("data/five_taxon_tree2.nwk", "data/five_taxon.fasta");
  auto [log_lik_trees_3_4, ignored_3_4] =
      ComputeExactMarginal("data/five_taxon_trees_3_4.nwk", "data/five_taxon.fasta");

  EigenVectorXd expected_log_lik_vector_at_rootsplits(3);
  expected_log_lik_vector_at_rootsplits << log_lik_tree_1, log_lik_tree_2,
      log_lik_trees_3_4;
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

TEST_CASE("GPInstance: inverted GPCSP probabilities") {
  // Note that just for fun, I have duplicated the first tree, but that doesn't matter
  // because we are looking at uniform over topological support.
  auto inst =
      GPInstanceOfFiles("data/five_taxon.fasta", "data/five_taxon_rooted_more_2.nwk");
  // See the DAG and the uniform probabilities at
  // https://github.com/phylovi/libsbn/issues/323#issuecomment-785410551
  const auto& dag = inst.GetDAG();
  EigenVectorXd normalized_sbn_parameters = dag.BuildUniformOnTopologicalSupportPrior();
  EigenVectorXd node_probabilities =
      dag.UnconditionalNodeProbabilities(normalized_sbn_parameters);
  EigenVectorXd correct_node_probabilities(15);
  correct_node_probabilities <<  //
      1.,                        // 0
      1.,                        // 1
      1.,                        // 2
      1.,                        // 3
      1.,                        // 4
      0.5,                       // 5
      0.25,                      // 6
      0.75,                      // 7
      0.25,                      // 8
      0.5,                       // 9
      0.25,                      // 10
      0.25,                      // 11
      0.5,                       // 12
      0.5,                       // 13
      0.25;                      // 14
  CheckVectorXdEquality(node_probabilities, correct_node_probabilities, 1e-12);

  EigenVectorXd inverted_probabilities =
      dag.InvertedGPCSPProbabilities(normalized_sbn_parameters, node_probabilities);
  EigenVectorXd correct_inverted_probabilities(24);
  correct_inverted_probabilities <<  //
                                     //
      1.,                            // 0 (rootsplit)
      1.,                            // 1 (rootsplit)
      1.,                            // 2 (rootsplit)
      1. / 3.,                       // 3
      0.5,                           // 4
      1.,                            // 5
      // We have the 0.5 coming from node 12, but that's split evenly between the two
      // descendants, so we have 0.25 from each. Thus even weights.
      0.5,      // 6
      1.,       // 7
      1.,       // 8
      1.,       // 9
      0.5,      // 10 (analogous to 6)
      2. / 3.,  // 11
      0.5,      // 12
      0.5,      // 13
      0.5,      // 14
      0.5,      // 15
      0.5,      // 16
      0.25,     // 17
      0.5,      // 18
      0.25,     // 19
      0.25,     // 20
      0.75,     // 21
      0.75,     // 22
      0.25;     // 23
  CheckVectorXdEquality(inverted_probabilities, correct_inverted_probabilities, 1e-12);
}

TEST_CASE("GPInstance: GenerateCompleteRootedTreeCollection") {
  const std::string fasta_path = "data/5-taxon-slice-of-ds1.fasta";
  auto inst =
      GPInstanceOfFiles(fasta_path, "data/5-taxon-only-rootward-uncertainty.nwk");
  EigenVectorXd branch_lengths(14);
  // The branch lengths contain the index of this GPCSP-indexed vector.
  branch_lengths << 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13.;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  // Because the branch lengths contain the GPCSP index, we can check that the indices
  // correspond to what we see in the GPCSP DAG in
  // https://github.com/phylovi/libsbn/issues/323#issuecomment-787451356
  CHECK_EQ(inst.GenerateCompleteRootedTreeCollection().Newick(),
           "((0:7,1:9):3,(2:11,(3:12,4:13):5):2):0;\n"
           "(1:10,(0:8,(2:11,(3:12,4:13):5):6):4):0;\n");
}

EigenVectorXd ClassicalLikelihoodOf(const std::string& tree_path,
                                    const std::string& fasta_path) {
  RootedSBNInstance sbn_instance("charlie");
  sbn_instance.ReadNewickFile(tree_path);
  sbn_instance.ProcessLoadedTrees();
  const Alignment alignment = Alignment::ReadFasta(fasta_path);
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  sbn_instance.SetAlignment(alignment);
  sbn_instance.PrepareForPhyloLikelihood(simple_specification, 1);

  std::vector<double> manual_log_likelihoods = sbn_instance.UnrootedLogLikelihoods();
  const double log_prior = log(1. / sbn_instance.tree_collection_.TreeCount());
  std::transform(manual_log_likelihoods.begin(), manual_log_likelihoods.end(),
                 manual_log_likelihoods.begin(),
                 [&log_prior](double log_like) { return log_like + log_prior; });
  return EigenVectorXdOfStdVectorDouble(manual_log_likelihoods);
}

// This is the simplest hybrid marginal that has tree uncertainty above and below the
// focal PCSP. Note that this test and the next one are set up so that the quartets
// reach far enough out that there is no uncertainty in the part of the tree outside of
// the quartet. In this case the hybrid marginal will be the same as the sum of
// classical likelihoods.
TEST_CASE("GPInstance: simplest hybrid marginal") {
  const std::string fasta_path = "data/7-taxon-slice-of-ds1.fasta";
  // See the DAG at
  // https://github.com/phylovi/libsbn/issues/323#issuecomment-783349464
  auto inst = GPInstanceOfFiles(fasta_path, "data/simplest-hybrid-marginal.nwk");
  auto& dag = inst.GetDAG();
  // Branch lengths generated from Python via
  // import random
  // [round(random.uniform(1e-6, 0.1), 3) for i in range(23)]
  EigenVectorXd branch_lengths(23);
  branch_lengths << 0.058, 0.044, 0.006, 0.099, 0.078, 0.036, 0.06, 0.073, 0.004, 0.041,
      0.088, 0.033, 0.043, 0.096, 0.027, 0.039, 0.043, 0.023, 0.064, 0.032, 0.03, 0.085,
      0.034;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  inst.PopulatePLVs();
  const std::string tree_path = "_ignore/simplest-hybrid-marginal-trees.nwk";
  inst.ExportAllGeneratedTrees(tree_path);

  // requests are printable to stdout if you're keen.
  auto request = dag.QuartetHybridRequestOf(12, false, 11);
  EigenVectorXd quartet_log_likelihoods =
      inst.GetEngine()->CalculateQuartetHybridLikelihoods(request);

  // Note that we aren't sorting likelihoods here, though we might have to do so for
  // more complex tests. I don't think that there's any guarantee that the hybrid log
  // likelihoods will be in the same order as the generated tree, but it worked here.
  EigenVectorXd manual_log_likelihoods = ClassicalLikelihoodOf(tree_path, fasta_path);
  CheckVectorXdEquality(quartet_log_likelihoods, manual_log_likelihoods, 1e-12);

  CHECK_EQ(request.IsFullyFormed(), true);
  CHECK_EQ(dag.QuartetHybridRequestOf(14, true, 13).IsFullyFormed(), false);
  CHECK_EQ(dag.QuartetHybridRequestOf(14, false, 0).IsFullyFormed(), false);
  CHECK_EQ(dag.QuartetHybridRequestOf(8, true, 4).IsFullyFormed(), false);
}

// This is a slightly more complex test, that has a rotation status of true, and has
// some paths through the DAG that aren't part of the hybrid marginal.
TEST_CASE("GPInstance: second simplest hybrid marginal") {
  const std::string fasta_path = "data/7-taxon-slice-of-ds1.fasta";
  // See the DAG at
  // https://github.com/phylovi/libsbn/issues/323#issuecomment-787963552
  auto inst = GPInstanceOfFiles(fasta_path, "data/second-simplest-hybrid-marginal.nwk");
  auto& dag = inst.GetDAG();
  // Branch lengths generated from Python via
  // import random
  // [round(random.uniform(1e-6, 0.1), 3) for i in range(32)]
  EigenVectorXd branch_lengths(32);
  branch_lengths << 0.09, 0.064, 0.073, 0.062, 0.051, 0.028, 0.077, 0.097, 0.089, 0.061,
      0.036, 0.049, 0.085, 0.01, 0.099, 0.027, 0.07, 0.023, 0.043, 0.056, 0.043, 0.026,
      0.058, 0.015, 0.093, 0.01, 0.011, 0.007, 0.022, 0.009, 0.037, 0.017;
  inst.GetEngine()->SetBranchLengths(branch_lengths);
  inst.PopulatePLVs();
  const std::string tree_path = "_ignore/simplest-hybrid-marginal-trees.nwk";
  inst.ExportAllGeneratedTrees(tree_path);

  auto request = dag.QuartetHybridRequestOf(12, true, 11);
  EigenVectorXd quartet_log_likelihoods =
      inst.GetEngine()->CalculateQuartetHybridLikelihoods(request);

  inst.LoadAllGeneratedTrees();
  // We restrict to only the trees that contain the DAG edge 6 (which goes between node
  // 12 and node 11). We get the bitset representation using inst.PrintGPCSPIndexer();
  inst.ExportTreesWithAPCSP("000000100111100001110", tree_path);
  EigenVectorXd manual_log_likelihoods = ClassicalLikelihoodOf(tree_path, fasta_path);
  CheckVectorXdEquality(quartet_log_likelihoods, manual_log_likelihoods, 1e-12);
}
