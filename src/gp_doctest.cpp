// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "include_doctest.hpp"

#include "sugar.hpp"
#include "combinatorics.hpp"
#include "gp_instance.hpp"
#include "phylo_model.hpp"
#include "reindexer.hpp"
#include "rooted_sbn_instance.hpp"
#include "stopwatch.hpp"
#include "tidy_subsplit_dag.hpp"
#include "nni_engine.hpp"
#include "pv_handler.hpp"
#include "topology_sampler.hpp"
#include "tp_engine.hpp"
#include "tp_choice_map.hpp"
#include "sankoff_matrix.hpp"
#include "sankoff_handler.hpp"
#include "dag_data.hpp"
#include "optimization.hpp"

using namespace GPOperations;  // NOLINT
using PLVType = PLVNodeHandler::PLVType;

// #350 remove all uses of GPCSP.

std::ostream& operator<<(std::ostream& os, EigenConstMatrixXdRef mx) {
  for (Eigen::Index i = 0; i < mx.rows(); i++) {
    for (Eigen::Index j = 0; j < mx.cols(); j++) {
      os << "[" << i << "," << j << "]: " << mx(i, j) << "\t";
    }
    os << std::endl;
  }
  return os;
}

// Let the "venus" node be the common ancestor of mars and saturn.
enum HelloGPCSP { jupiter, mars, saturn, venus, rootsplit, root };

// *** GPInstances used for testing ***

GPInstance GPInstanceOfFiles(
    const std::string& fasta_path, const std::string& newick_path,
    const std::string mmap_filepath = std::string("_ignore/mmapped_pv.data"),
    const bool use_gradients = false) {
  GPInstance inst(mmap_filepath);
  inst.ReadFastaFile(fasta_path);
  inst.ReadNewickFile(newick_path, false);
  inst.MakeDAG();
  inst.MakeGPEngine();
  inst.UseGradientOptimization(use_gradients);
  return inst;
}

// Our tree is (see check below)
// (jupiter:0.113,(mars:0.15,saturn:0.1)venus:0.22):0.;
// You can see a helpful diagram at
// https://github.com/phylovi/bito/issues/349#issuecomment-898672399
GPInstance MakeHelloGPInstance(const std::string& fasta_path) {
  auto inst = GPInstanceOfFiles(fasta_path, "data/hello_rooted.nwk");
  EigenVectorXd branch_lengths(5);
  // Order set by HelloGPCSP.
  branch_lengths << 0, 0.22, 0.113, 0.15, 0.1;
  inst.GetGPEngine().SetBranchLengths(branch_lengths);
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
// https://github.com/phylovi/bito/issues/391#issuecomment-1169048568
GPInstance MakeDS1Reduced5Instance() {
  auto inst = GPInstanceOfFiles("data/ds1-reduced-5.fasta", "data/ds1-reduced-5.nwk");
  return inst;
}

GPInstance MakeFluAGPInstance(double rescaling_threshold) {
  auto inst = GPInstanceOfFiles("data/fluA.fa", "data/fluA.tree");
  inst.MakeGPEngine(rescaling_threshold);
  inst.GetGPEngine().SetBranchLengthsToConstant(0.01);
  return inst;
}

TEST_CASE("DAGSummaryStatistics") {
  auto inst = MakeHelloGPInstanceTwoTrees();
  StringSizeMap summaries = {{"edge_count", 10}, {"node_count", 8}};
  CHECK(summaries == inst.DAGSummaryStatistics());
}

EigenVectorXd MakeHelloGPInstanceMarginalLikelihoodTestBranchLengths() {
  EigenVectorXd hello_gp_optimal_branch_lengths(10);
  hello_gp_optimal_branch_lengths << 1, 1, 0.066509261, 0.00119570257, 0.00326456973,
      0.0671995398, 0.203893516, 0.204056242, 0.0669969961, 0.068359082;

  return hello_gp_optimal_branch_lengths;
}

TEST_CASE("GPInstance: straightforward classical likelihood calculation") {
  auto inst = MakeHelloGPInstance();
  auto& engine = inst.GetGPEngine();

  inst.PopulatePLVs();
  inst.ComputeLikelihoods();

  EigenVectorXd realized_log_likelihoods =
      inst.GetGPEngine().GetPerGPCSPLogLikelihoods();
  CheckVectorXdEquality(-84.77961943, realized_log_likelihoods, 1e-6);

  CHECK_LT(fabs(engine.GetLogMarginalLikelihood() - -84.77961943), 1e-6);
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
  sbn_instance.ReadNewickFile(newick_path, false);
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
      Assert(Bitset(gp_string.substr(0, gp_string.find('|') - 1)).None() ||
                 Bitset(gp_string.substr(gp_string.rfind('|') + 1)).None(),
             "Missing an internal node in CheckExactMapVsGPVector.");
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
  inst.EstimateBranchLengths(0.00001, 100, true);

  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  inst.ComputeMarginalLikelihood();
  std::string tree_path = "_ignore/test_marginal_trees.nwk";
  const auto trees = inst.CurrentlyLoadedTreesWithGPBranchLengths();
  trees.ToNewickFile(tree_path);

  auto [exact_log_likelihood, exact_per_pcsp_log_marginal] =
      ComputeExactMarginal(tree_path, fasta_path);
  double gp_marginal_log_likelihood = inst.GetGPEngine().GetLogMarginalLikelihood();
  auto gp_per_pcsp_log_marginal =
      inst.PrettyIndexedPerGPCSPComponentsOfFullLogMarginal();

  double tolerance = 1e-6;
  if (fabs(gp_marginal_log_likelihood - exact_log_likelihood) > tolerance) {
    std::cout << "gp_marginal_log_likelihood: " << gp_marginal_log_likelihood
              << std::endl;
    std::cout << "exact_log_likelihood: " << exact_log_likelihood << std::endl;
  }
  CHECK_LT(fabs(gp_marginal_log_likelihood - exact_log_likelihood), tolerance);
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
  // https://github.com/phylovi/bito/issues/391#issuecomment-1169053191
  TestCompositeMarginal(
      GPInstanceOfFiles(fasta_path, "data/simplest-hybrid-marginal-all-trees.nwk"),
      fasta_path);
}

TEST_CASE("GPInstance: gradient calculation") {
  auto inst = MakeHelloGPInstanceSingleNucleotide();
  auto& engine = inst.GetGPEngine();

  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  NodeId rootsplit_id = rootsplit;
  NodeId child_id = jupiter;
  NodeId rootsplit_jupiter_idx = 2;
  size_t hello_node_count_without_dag_root_node = 5;

  PVId leafward_idx = PLVNodeHandler::GetPVIndex(
      PLVType::P, child_id, hello_node_count_without_dag_root_node);
  PVId rootward_idx = PLVNodeHandler::GetPVIndex(
      PLVType::RLeft, rootsplit_id, hello_node_count_without_dag_root_node);
  OptimizeBranchLength op{leafward_idx.value_, rootward_idx.value_,
                          rootsplit_jupiter_idx.value_};
  DoublePair log_lik_and_derivative = engine.LogLikelihoodAndDerivative(op);
  // Expect log lik: -4.806671945.
  // Expect log lik derivative: -0.6109379521.
  CHECK_LT(fabs(log_lik_and_derivative.first - -4.806671945), 1e-6);
  CHECK_LT(fabs(log_lik_and_derivative.second - -0.6109379521), 1e-6);
}

TEST_CASE("GPInstance: multi-site gradient calculation") {
  auto inst = MakeHelloGPInstance();
  auto& engine = inst.GetGPEngine();

  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  NodeId rootsplit_id = rootsplit;
  NodeId child_id = jupiter;
  NodeId rootsplit_jupiter_idx = 2;
  size_t hello_node_count_without_dag_root_node = 5;

  PVId leafward_idx = PLVNodeHandler::GetPVIndex(
      PLVType::P, child_id, hello_node_count_without_dag_root_node);
  PVId rootward_idx = PLVNodeHandler::GetPVIndex(
      PLVType::RLeft, rootsplit_id, hello_node_count_without_dag_root_node);
  OptimizeBranchLength op{leafward_idx.value_, rootward_idx.value_,
                          rootsplit_jupiter_idx.value_};
  std::tuple<double, double, double> log_lik_and_derivatives =
      engine.LogLikelihoodAndFirstTwoDerivatives(op);
  // Expect log likelihood: -84.77961943.
  // Expect log llh first derivative: -18.22479569.
  // Expect log llh second derivative: -5.4460787413.
  CHECK_LT(fabs(std::get<0>(log_lik_and_derivatives) - -84.77961943), 1e-6);
  CHECK_LT(fabs(std::get<1>(log_lik_and_derivatives) - -18.22479569), 1e-6);
  CHECK_LT(fabs(std::get<2>(log_lik_and_derivatives) - -5.4460787413), 1e-6);
}

// We are outputting the branch length for PCSP 100-011-001
// which has a true branch length of 0.0694244266
double ObtainBranchLengthWithOptimization(const OptimizationMethod method,
                                          bool is_quiet = true) {
  GPInstance inst = MakeHelloGPInstance();
  GPEngine& engine = inst.GetGPEngine();
  engine.SetOptimizationMethod(method);
  inst.EstimateBranchLengths(0.0001, 100, is_quiet);

  inst.MakeDAG();
  GPDAG& dag = inst.GetDAG();
  EdgeId default_index = EdgeId(dag.EdgeCountWithLeafSubsplits());
  Bitset gpcsp_bitset = Bitset("100011001");
  EdgeId index =
      AtWithDefault(dag.BuildEdgeIndexer(), gpcsp_bitset, default_index.value_);
  return inst.GetGPEngine().GetBranchLengths()(index.value_);
}

TEST_CASE("GPInstance: Gradient-based optimization with Newton's Method") {
  double nongradient_length =
      ObtainBranchLengthWithOptimization(OptimizationMethod::BrentOptimization);
  double gradient_length =
      ObtainBranchLengthWithOptimization(OptimizationMethod::NewtonOptimization);

  double true_length = 0.0694244266;
  double nongrad_diff = abs(nongradient_length - true_length);
  double grad_diff = abs(gradient_length - true_length);
  double tol = 1e-6;

  if (grad_diff > nongrad_diff || grad_diff > tol) {
    std::cout << "nongrad_diff: " << nongrad_diff << std::endl;
    std::cout << "grad_diff: " << grad_diff << std::endl;
    std::cout << "nongradient_length: " << nongradient_length << std::endl;
    std::cout << "gradient_length: " << gradient_length << std::endl;
    std::cout << "true_length: " << true_length << std::endl;
  }
  CHECK_LT(grad_diff, nongrad_diff);
  CHECK_LT(grad_diff, tol);
}

double MakeAndRunFluAGPInstance(double rescaling_threshold) {
  auto inst = MakeFluAGPInstance(rescaling_threshold);
  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  return inst.GetGPEngine().GetLogMarginalLikelihood();
}

TEST_CASE("GPInstance: rescaling") {
  double difference = MakeAndRunFluAGPInstance(GPEngine::default_rescaling_threshold_) -
                      MakeAndRunFluAGPInstance(1e-4);
  CHECK_LT(fabs(difference), 1e-10);
}

StringDoubleMap StringDoubleMapOfStringDoubleVector(StringDoubleVector vect) {
  StringDoubleMap m;
  for (const auto& [str, x] : vect) {
    SafeInsert(m, str, x);
  }
  return m;
}

TEST_CASE("GPInstance: gather and hotstart branch lengths") {
  // » nw_topology data/hotstart_bootstrap_sample.nwk | nw_order - | sort | uniq -c
  // 1 (outgroup,(((z0,z1),z2),z3));
  // 33 (outgroup,((z0,z1),(z2,z3)));
  const std::string tree_path = "data/hotstart_bootstrap_sample.nwk";
  GPInstance inst("_ignore/mmapped_pv.data");
  // This is just a dummy fasta file, which is required to make an Engine.
  inst.ReadFastaFile("data/hotstart.fasta");
  inst.ReadNewickFile(tree_path, false);
  inst.MakeGPEngine();

  // We are going to verify correct assignment of the PCSP with sister z2, z3 and
  // children z0, z1, which only appears in the tree (outgroup,((z0,z1),(z2,z3))).
  // Vector of taxon names: [outgroup, z2, z3, z1, z0]
  // So, this below is the desired GPCSP: 01100|00011|00001

  // These branch lengths are obtained by excluding (outgroup,(((z0,z1),z2),z3)) (which
  // doesn't have this PCSP) and grabbing the rest of the branch lengths.
  //
  // We will first test that gather branch lengths is collecting the correct set of
  // branches, and then we will test whether hot start is accurately calculating the
  // mean of these branches.

  EigenVectorXd expected_bls_internal(33);
  expected_bls_internal << 0.1175370000, 0.1175750000, 0.1195780000, 0.0918962000,
      0.0918931000, 0.1192590000, 0.0906988000, 0.0906972000, 0.0905154000,
      0.0903663000, 0.1245620000, 0.1244890000, 0.1245050000, 0.1245550000,
      0.1245680000, 0.1248920000, 0.1248490000, 0.1164070000, 0.1164110000,
      0.1164120000, 0.1245670000, 0.1245650000, 0.1245670000, 0.1245670000,
      0.1240790000, 0.1242540000, 0.1242160000, 0.1242560000, 0.1892030000,
      0.1894900000, 0.1895430000, 0.1896900000, 0.1905710000;

  SizeDoubleVectorMap branch_lengths_from_sample = inst.GatherBranchLengths();
  EigenVectorXd gathered_bls =
      EigenVectorXdOfStdVectorDouble(branch_lengths_from_sample[4]);
  CheckVectorXdEquality(expected_bls_internal, gathered_bls, 1e-6);

  double true_mean_internal = expected_bls_internal.array().mean();
  inst.HotStartBranchLengths();
  CHECK_LT(fabs(true_mean_internal - inst.GetGPEngine().GetBranchLengths()(4)), 1e-8);
  // We also want to verify correct assignment for a pendant branch length.
  // Specifically, we are looking at the pendant branch length for z2 with sister z3. So
  // the desired GPCSP is 0010001000|0000000000. This corresponds to branch length index
  // 8, and is also found by excluding (outgroup, (((z0,z1),z2),z3)), which does not
  // have this PCSP.
  EigenVectorXd expected_bls_pendant(33);
  expected_bls_pendant << 0.0903520000, 0.0903100000, 0.0911710000, 0.0906700000,
      0.0906680000, 0.0907450000, 0.0884430000, 0.0883790000, 0.0909010000,
      0.0865700000, 0.0999870000, 0.0999920000, 0.0999680000, 0.0999430000,
      0.0999610000, 0.0902300000, 0.0902700000, 0.0905340000, 0.0908440000,
      0.0901110000, 0.0898580000, 0.0898570000, 0.0909610000, 0.0898660000,
      0.0906510000, 0.0906750000, 0.0906480000, 0.0906100000, 0.0894660000,
      0.0904620000, 0.0893220000, 0.0902220000, 0.0902000000;
  double true_mean_pendant = expected_bls_pendant.array().mean();
  CHECK_LT(fabs(true_mean_pendant - inst.GetGPEngine().GetBranchLengths()(8)), 1e-8);
}

TEST_CASE("GPInstance: take first branch length") {
  const std::string tree_path = "data/hotstart_bootstrap_sample.nwk";
  GPInstance inst("_ignore/mmapped_pv.data");
  // This is just a dummy fasta file, which is required to make an Engine.
  inst.ReadFastaFile("data/hotstart.fasta");
  inst.ReadNewickFile(tree_path, false);
  inst.MakeGPEngine();
  inst.TakeFirstBranchLength();
  auto branch_length_map =
      StringDoubleMapOfStringDoubleVector(inst.PrettyIndexedBranchLengths());

  // Check at the internal branches
  CHECK_EQ(0.191715, branch_length_map.at("01000|00011|00001"));     // pcsp index 1
  CHECK_EQ(0.117537, branch_length_map.at("01100|00011|00001"));     // pcsp index 2
  CHECK_EQ(0.0874183, branch_length_map.at("00011|01100|00100"));    // pcsp index 3
  CHECK_EQ(0.129921, branch_length_map.at("10000|01111|00011"));     // pcsp index 4
  CHECK_EQ(0.15936, branch_length_map.at("10000|01111|00100"));      // pcsp index 5
  CHECK_EQ(0.000813992, branch_length_map.at("00100|01011|00011"));  // pcsp index 6

  // Check at the branches ending in leaves
  CHECK_EQ(0.129921, branch_length_map.at("01111|10000|00000"));  // pcsp index 7
  CHECK_EQ(0.090352, branch_length_map.at("00100|01000|00000"));  // pcsp index 8
  CHECK_EQ(0.099922, branch_length_map.at("00011|01000|00000"));  // pcsp index 9
  CHECK_EQ(0.112125, branch_length_map.at("01000|00100|00000"));  // pcsp index 10
  CHECK_EQ(0.104088, branch_length_map.at("01011|00100|00000"));  // pcsp index 11
  CHECK_EQ(0.113775, branch_length_map.at("00001|00010|00000"));  // pcsp index 12
  CHECK_EQ(0.081634, branch_length_map.at("00010|00001|00000"));  // pcsp index 13
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
  size_t length = inst.GetGPEngine().GetLogLikelihoodMatrix().rows();
  const EigenVectorXd log_likelihoods1 =
      inst.GetGPEngine().GetPerGPCSPLogLikelihoods(0, length);
  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  const EigenVectorXd log_likelihoods2 = inst.GetGPEngine().GetPerGPCSPLogLikelihoods();
  CheckVectorXdEquality(log_likelihoods1, log_likelihoods2, 1e-6);
}

TEST_CASE("GPInstance: SBN root split probabilities on five taxa") {
  auto inst = MakeFiveTaxonInstance();
  inst.GetGPEngine().SetBranchLengthsToConstant(0.1);
  inst.PopulatePLVs();
  // We need to call ComputeLikelihoods to populate the likelihood matrix.
  // Note: EstimateBranchLengths doesn't populate the likelihood matrix.
  inst.ComputeLikelihoods();

  EigenVectorXd log_likelihood_vector = inst.GetGPEngine().GetPerGPCSPLogLikelihoods();

  // Let s be a subsplit and k be the site. Then,
  // log_likelihood_matrix.row(s)[k] =
  //    \log \sum_{\tau : s \in \tau} q(\tau) P(y_k | \tau),
  // log_likelihood_vector[s] =
  //    \sum_{k=1}^{K} \log \sum_{\tau : s \in \tau} q(\tau) P(y_k | \tau).
  // To test this, we are going to compute P(y_k | \tau) for {\tau : s \in \tau} and
  // multiply this by q(\tau) = 1/4 since we are assuming a uniform prior.

  // The collection of trees that we are looking at has 3 rootplits where one root
  // split generates two trees and the other 2 root splits generating one tree each
  // for the total of 4 trees.

  // We will compare the values against the 3 rootsplits, since we cannot assume
  // the ordering due to different implementation of the map, we will sort the values
  // before comparison.

  auto [log_lik_tree_1, ignored_1] =
      ComputeExactMarginal("data/five_taxon_tree1.nwk", "data/five_taxon.fasta");
  std::ignore = ignored_1;
  auto [log_lik_tree_2, ignored_2] =
      ComputeExactMarginal("data/five_taxon_tree2.nwk", "data/five_taxon.fasta");
  std::ignore = ignored_2;
  auto [log_lik_trees_3_4, ignored_3_4] =
      ComputeExactMarginal("data/five_taxon_trees_3_4.nwk", "data/five_taxon.fasta");
  std::ignore = ignored_3_4;

  EigenVectorXd expected_log_lik_vector_at_rootsplits(3);
  expected_log_lik_vector_at_rootsplits << log_lik_tree_1, log_lik_tree_2,
      log_lik_trees_3_4;
  EigenVectorXd realized_log_lik_vector_at_rootsplits =
      log_likelihood_vector.segment(0, 3);
  CheckVectorXdEqualityAfterSorting(realized_log_lik_vector_at_rootsplits,
                                    expected_log_lik_vector_at_rootsplits, 1e-6);

  inst.EstimateSBNParameters();
  EigenVectorXd realized_q = inst.GetGPEngine().GetSBNParameters().segment(0, 3);
  // The expected values for the SBN parameters: q[s] \propto log_lik[s] +
  // log_prior[s]. The SBN params are initialized so that we get a uniform
  // distribution over the trees. For the rootsplits, the values are (1/4, 1/4, 2/4)
  // corresponding to the entries in expected_log_lik_vector_at_rootsplits.
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
  inst.GetGPEngine().SetBranchLengths(branch_lengths);
  auto trees = inst.CurrentlyLoadedTreesWithGPBranchLengths();
  CHECK_EQ(trees.Newick(), "(jupiter:0.2,(mars:0.3,saturn:0.4):0.1):0;\n");
}

TEST_CASE("GPInstance: CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths") {
  GPInstance inst("_ignore/mmapped_pv.data");
  inst.ReadFastaFile("data/five_taxon.fasta");
  inst.ReadNewickFile("data/five_taxon_rooted_more.nwk", false);
  inst.MakeGPEngine();
  inst.GetGPEngine().SetBranchLengthsToConstant(0.9);
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
  // 0000|1111|0001,    0
  // 0000|1111|0011,    1
  // 0001|1110|0110,    2
  // 0001|1110|0010,    3
  auto support = inst.GetDAG().BuildUniformOnTopologicalSupportPrior();
  CHECK_LT(fabs(support[0] - 2. / 3.), 1e-10);
  CHECK_LT(fabs(support[1] - 1. / 3.), 1e-10);
  CHECK_LT(fabs(support[2] - 1. / 2.), 1e-10);
  CHECK_LT(fabs(support[3] - 1. / 2.), 1e-10);
  auto all = inst.GetDAG().BuildUniformOnAllTopologiesPrior();
  // There are 15 topologies on 4 taxa.
  // There are 3 topologies on 3 taxa, so there are 3 topologies with rootsplit
  // 0001|1110.
  CHECK_LT(fabs(all[0] - 3. / 15.), 1e-10);
  // There is only 1 topology with rootsplit 0011|1100.
  CHECK_LT(fabs(all[1] - 1. / 15.), 1e-10);
  // There are 3 topologies on 3 taxa.
  CHECK_LT(fabs(all[2] - 1. / 3.), 1e-10);
  CHECK_LT(fabs(all[3] - 1. / 3.), 1e-10);
}

TEST_CASE("GPInstance: inverted GPCSP probabilities") {
  // Note that just for fun, I have duplicated the first tree, but that doesn't matter
  // because we are looking at uniform over topological support.
  auto inst =
      GPInstanceOfFiles("data/five_taxon.fasta", "data/five_taxon_rooted_more_2.nwk");
  // See the DAG and the uniform probabilities at
  // https://github.com/phylovi/bito/issues/391#issuecomment-1168046752
  const auto& dag = inst.GetDAG();
  EigenVectorXd normalized_sbn_parameters = dag.BuildUniformOnTopologicalSupportPrior();
  EigenVectorXd node_probabilities =
      dag.UnconditionalNodeProbabilities(normalized_sbn_parameters);
  EigenVectorXd correct_node_probabilities(16);
  correct_node_probabilities <<  //
      1.,                        // 0
      1.,                        // 1
      1.,                        // 2
      1.,                        // 3
      1.,                        // 4
      0.75,                      // 5
      0.5,                       // 6
      0.25,                      // 7
      0.25,                      // 8
      0.5,                       // 9
      0.25,                      // 10
      0.25,                      // 11
      0.5,                       // 12
      0.5,                       // 13
      0.25,                      // 14
      1.;                        // 15 (DAG root node)
  CheckVectorXdEquality(node_probabilities, correct_node_probabilities, 1e-12);

  EigenVectorXd inverted_probabilities =
      dag.InvertedGPCSPProbabilities(normalized_sbn_parameters, node_probabilities);
  EigenVectorXd correct_inverted_probabilities(24);
  correct_inverted_probabilities <<  //
                                     //
      1.,                            // 0 (rootsplit)
      1.,                            // 1 (rootsplit)
      1.,                            // 2 (rootsplit)
      1.,                            // 3
      1.,                            // 4
      2. / 3.,                       // 5
      0.5,                           // 6
      0.5,                           // 7
      // We have 0.5 from node 9, but that's split proportionally to the probability
      // of each potential parent. Nodes 12 and 14 are equally likely parents of node 9,
      // so we have 0.5 for the inverted PCSP probability.
      0.5,      // 8
      1.,       // 9
      1.,       // 10
      0.5,      // 11 (analogous to 8)
      1. / 3.,  // 12
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
  inst.GetGPEngine().SetBranchLengths(branch_lengths);
  // Because the branch lengths contain the GPCSP index, we can check that the indices
  // correspond to what we see in the GPCSP DAG in
  // https://github.com/phylovi/bito/issues/391#issuecomment-1168048090
  CHECK_EQ(inst.GenerateCompleteRootedTreeCollection().Newick(),
           "((0:7,1:9):3,(2:11,(3:12,4:13):2):6):0;\n"
           "(1:10,(0:8,(2:11,(3:12,4:13):2):5):4):0;\n");
}

EigenVectorXd ClassicalLikelihoodOf(const std::string& tree_path,
                                    const std::string& fasta_path) {
  RootedSBNInstance sbn_instance("charlie");
  sbn_instance.ReadNewickFile(tree_path, false);
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
// reach far enough out that there is no uncertainty in the part of the tree outside
// of the quartet. In this case the hybrid marginal will be the same as the sum of
// classical likelihoods.
TEST_CASE("GPInstance: simplest hybrid marginal") {
  const std::string fasta_path = "data/7-taxon-slice-of-ds1.fasta";
  // See the DAG at
  // https://github.com/phylovi/bito/issues/391#issuecomment-1169053191
  auto inst = GPInstanceOfFiles(fasta_path, "data/simplest-hybrid-marginal.nwk");
  auto& dag = inst.GetDAG();
  // Branch lengths generated from Python via
  // import random
  // [round(random.uniform(1e-6, 0.1), 3) for i in range(23)]
  EigenVectorXd branch_lengths(23);
  branch_lengths << 0.058, 0.044, 0.006, 0.099, 0.078, 0.036, 0.06, 0.073, 0.004, 0.041,
      0.088, 0.033, 0.043, 0.096, 0.027, 0.039, 0.043, 0.023, 0.064, 0.032, 0.03, 0.085,
      0.034;
  inst.GetGPEngine().SetBranchLengths(branch_lengths);
  inst.PopulatePLVs();
  const std::string tree_path = "_ignore/simplest-hybrid-marginal-trees.nwk";
  inst.ExportAllGeneratedTrees(tree_path);

  // requests are printable to stdout if you're keen.
  auto request = dag.QuartetHybridRequestOf(NodeId(12), false, NodeId(11));
  EigenVectorXd quartet_log_likelihoods =
      inst.GetGPEngine().CalculateQuartetHybridLikelihoods(request);

  // Note that we aren't sorting likelihoods here, though we might have to do so for
  // more complex tests. I don't think that there's any guarantee that the hybrid log
  // likelihoods will be in the same order as the generated tree, but it worked here.
  EigenVectorXd manual_log_likelihoods = ClassicalLikelihoodOf(tree_path, fasta_path);
  CheckVectorXdEquality(quartet_log_likelihoods, manual_log_likelihoods, 1e-12);

  CHECK_EQ(request.IsFullyFormed(), true);
  CHECK_EQ(dag.QuartetHybridRequestOf(NodeId(14), true, NodeId(13)).IsFullyFormed(),
           false);
  CHECK_EQ(dag.QuartetHybridRequestOf(NodeId(14), false, NodeId(0)).IsFullyFormed(),
           false);
  CHECK_EQ(dag.QuartetHybridRequestOf(NodeId(8), true, NodeId(4)).IsFullyFormed(),
           false);
}

// This is a slightly more complex test, that has a rotation status of true, and has
// some paths through the DAG that aren't part of the hybrid marginal.
TEST_CASE("GPInstance: second simplest hybrid marginal") {
  const std::string fasta_path = "data/7-taxon-slice-of-ds1.fasta";
  // See the DAG at
  // https://github.com/phylovi/bito/issues/391#issuecomment-1169056581
  auto inst = GPInstanceOfFiles(fasta_path, "data/second-simplest-hybrid-marginal.nwk");
  auto& dag = inst.GetDAG();
  // Branch lengths generated from Python via
  // import random
  // [round(random.uniform(1e-6, 0.1), 3) for i in range(32)]
  EigenVectorXd branch_lengths(32);
  branch_lengths << 0.09, 0.064, 0.073, 0.062, 0.051, 0.028, 0.077, 0.097, 0.089, 0.061,
      0.036, 0.049, 0.085, 0.01, 0.099, 0.027, 0.07, 0.023, 0.043, 0.056, 0.043, 0.026,
      0.058, 0.015, 0.093, 0.01, 0.011, 0.007, 0.022, 0.009, 0.037, 0.017;
  inst.GetGPEngine().SetBranchLengths(branch_lengths);
  inst.PopulatePLVs();
  const std::string tree_path = "_ignore/simplest-hybrid-marginal-trees.nwk";
  inst.ExportAllGeneratedTrees(tree_path);

  auto edge = dag.GetDAGEdge(EdgeId(2));
  auto request = dag.QuartetHybridRequestOf(NodeId(edge.GetParent()), true,
                                            NodeId(edge.GetChild()));
  EigenVectorXd quartet_log_likelihoods =
      inst.GetGPEngine().CalculateQuartetHybridLikelihoods(request);

  inst.LoadAllGeneratedTrees();
  // We restrict to only the trees that contain the DAG edge 6 (which goes between
  // node 12 and node 11). We get the bitset representation using
  // inst.PrintGPCSPIndexer();
  inst.ExportTreesWithAPCSP("000000100111100001110", tree_path);
  EigenVectorXd manual_log_likelihoods = ClassicalLikelihoodOf(tree_path, fasta_path);
  CheckVectorXdEquality(quartet_log_likelihoods, manual_log_likelihoods, 1e-12);
}

TEST_CASE("GPInstance: test GPCSP indexes") {
  const std::string fasta_path = "data/7-taxon-slice-of-ds1.fasta";
  auto inst = GPInstanceOfFiles(fasta_path, "data/simplest-hybrid-marginal.nwk");
  auto& dag = inst.GetDAG();
  dag.TopologicalEdgeTraversal([&dag](NodeId parent_id, bool is_edge_on_left,
                                      NodeId child_id, EdgeId gpcsp_idx) {
    CHECK_EQ(dag.GetEdgeIdx(parent_id, child_id), gpcsp_idx);
  });
}

// ** SubsplitDAG tests **

template <typename T>
std::vector<T> ConvertIdVector(const SizeVector& vec_in) {
  std::vector<T> vec_out;
  for (const auto i : vec_in) {
    vec_out.push_back(T(i));
  }
  return vec_out;
}

TEST_CASE("SubsplitDAG: test rootsplits") {
  const std::string fasta_path = "data/7-taxon-slice-of-ds1.fasta";
  auto inst = GPInstanceOfFiles(fasta_path, "data/simplest-hybrid-marginal.nwk");
  inst.SubsplitDAGToDot("_ignore/outtest.dot", true);
  auto& dag = inst.GetDAG();
  for (const auto& rootsplit_id : dag.GetRootsplitNodeIds()) {
    const auto rootsplit_node = dag.GetDAGNode(NodeId(rootsplit_id));
    CHECK(rootsplit_node.IsRootsplit());
  }
}

// See diagram at https://github.com/phylovi/bito/issues/391#issuecomment-1168046752.
TEST_CASE("SubsplitDAG: IsValidAddNodePair tests") {
  const std::string fasta_path = "data/five_taxon.fasta";
  auto inst = GPInstanceOfFiles(fasta_path, "data/five_taxon_rooted_more_2.nwk");
  auto& dag = inst.GetDAG();

  // Nodes are not adjacent (12|34 and 2|4).
  CHECK_FALSE(dag.IsValidAddNodePair(Bitset::Subsplit("01100", "00011"),
                                     Bitset::Subsplit("00100", "00001")));
  // Nodes have 6 taxa while the DAG has 5 (12|34 and 1|2).
  CHECK_FALSE(dag.IsValidAddNodePair(Bitset::Subsplit("011000", "000110"),
                                     Bitset::Subsplit("010000", "001000")));
  // Parent node does not have a parent (12|3 and 1|2).
  CHECK_FALSE(dag.IsValidAddNodePair(Bitset::Subsplit("01100", "00010"),
                                     Bitset::Subsplit("01000", "00100")));
  // Left clade of the parent node does not have a child (02|134 and
  // 1|34).
  CHECK_FALSE(dag.IsValidAddNodePair(Bitset::Subsplit("10100", "01011"),
                                     Bitset::Subsplit("01000", "00011")));
  // Left clade of the child node does not have a child (0123|4 and
  // 023|1).
  CHECK_FALSE(dag.IsValidAddNodePair(Bitset::Subsplit("11110", "00001"),
                                     Bitset::Subsplit("10110", "01000")));
  // Right clade of the child node does not have a child (0123|4 and
  // 0|123).
  CHECK_FALSE(dag.IsValidAddNodePair(Bitset::Subsplit("11110", "00001"),
                                     Bitset::Subsplit("10000", "01110")));
  // Valid new node pair (0123|4 and 012|3).
  CHECK(dag.IsValidAddNodePair(Bitset::Subsplit("11110", "00001"),
                               Bitset::Subsplit("11100", "00010")));
}

// See diagram at https://github.com/phylovi/bito/issues/391#issuecomment-1168059272.
TEST_CASE("SubsplitDAG: AddNodePair tests") {
  const std::string fasta_path = "data/five_taxon.fasta";
  auto inst = GPInstanceOfFiles(fasta_path, "data/five_taxon_rooted_more_2.nwk");
  auto& dag = inst.GetDAG();

  // Check that AddNodePair throws if node pair is invalid (12|34 and 2|4).
  CHECK_THROWS(dag.AddNodePair(Bitset::Subsplit("01100", "00011"),
                               Bitset::Subsplit("00100", "00001")));
  // Add 2|34 and 3|4, which are both already in the DAG.
  // Check that AddNodePair returns empty added_node_ids and added_edge_idxs
  // and that node_reindexer and edge_reindexer are the identity reindexers.
  auto node_addition_result = dag.AddNodePair(Bitset::Subsplit("00100", "00011"),
                                              Bitset::Subsplit("00010", "00001"));
  CHECK(node_addition_result.added_node_ids.empty());
  CHECK(node_addition_result.added_edge_idxs.empty());
  CHECK_EQ(node_addition_result.node_reindexer, Reindexer::IdentityReindexer(16));
  CHECK_EQ(node_addition_result.edge_reindexer, Reindexer::IdentityReindexer(24));
  // Before adding any nodes.
  size_t prev_node_count = dag.NodeCount();
  size_t prev_edge_count = dag.EdgeCountWithLeafSubsplits();
  size_t prev_topology_count = dag.TopologyCount();
  // Add nodes 24|3 and 2|4.
  Bitset parent_subsplit = Bitset::Subsplit("00101", "00010");
  Bitset child_subsplit = Bitset::Subsplit("00100", "00001");
  node_addition_result = dag.AddNodePair(parent_subsplit, child_subsplit);
  // Check that the node count and edge count was updated.
  CHECK_EQ(dag.NodeCount(), prev_node_count + 2);
  CHECK_EQ(dag.EdgeCountWithLeafSubsplits(), prev_edge_count + 6);
  // Check that both nodes now exist.
  CHECK(dag.ContainsNode(parent_subsplit));
  CHECK(dag.ContainsNode(child_subsplit));
  // Check that all necessary edges were created.
  const auto parent_node = dag.GetDAGNode(dag.GetDAGNodeId(parent_subsplit));
  const auto child_node = dag.GetDAGNode(dag.GetDAGNodeId(child_subsplit));
  std::map<bool, SizeVector> correct_parents_of_parent{{true, {}}, {false, {14, 16}}};
  std::map<bool, SizeVector> parents_of_parent{{true, parent_node.GetLeftRootward()},
                                               {false, parent_node.GetRightRootward()}};
  CHECK_EQ(parents_of_parent, correct_parents_of_parent);
  std::map<bool, SizeVector> children_of_parent{
      {true, parent_node.GetLeftLeafward()}, {false, parent_node.GetRightLeafward()}};
  std::map<bool, SizeVector> correct_children_of_parent{{true, {12}}, {false, {3}}};
  CHECK_EQ(children_of_parent, correct_children_of_parent);
  std::map<bool, SizeVector> parents_of_children{
      {true, child_node.GetLeftRootward()}, {false, child_node.GetRightRootward()}};
  std::map<bool, SizeVector> correct_parents_of_children{{true, {13}}, {false, {}}};
  CHECK_EQ(parents_of_children, correct_parents_of_children);
  std::map<bool, SizeVector> children_of_child{{true, child_node.GetLeftLeafward()},
                                               {false, child_node.GetRightLeafward()}};
  std::map<bool, SizeVector> correct_children_of_child{{true, {2}}, {false, {4}}};
  CHECK_EQ(children_of_child, correct_children_of_child);
  // Check that node_reindexer and edge_reindexer are correct.
  Reindexer correct_node_reindexer(
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 12, 13});
  CHECK_EQ(node_addition_result.node_reindexer, correct_node_reindexer);
  Reindexer correct_edge_reindexer({0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                    11, 12, 14, 15, 16, 17, 18, 19, 20, 21,
                                    22, 23, 24, 25, 26, 27, 28, 29, 13, 10});
  CHECK_EQ(node_addition_result.edge_reindexer, correct_edge_reindexer);
  // Check that added_node_ids and added_edge_idxs are correct.
  NodeIdVector correct_added_node_ids = ConvertIdVector<NodeId>({12, 13});
  CHECK_EQ(node_addition_result.added_node_ids, correct_added_node_ids);
  EdgeIdVector correct_added_edge_idxs =
      ConvertIdVector<EdgeId>({26, 27, 28, 29, 13, 10});
  CHECK_EQ(node_addition_result.added_edge_idxs, correct_added_edge_idxs);
  // Check that `dag_nodes` was updated (node 12 -> 14).
  const auto& node_14 = dag.GetDAGNode(NodeId(14));
  CHECK_EQ(node_14.GetBitset().ToString(), "0100000111");
  // Check that node fields were updated correctly.
  const auto& right_parents_14 = node_14.GetRightRootward();
  const auto& right_children_14 = node_14.GetRightLeafward();
  CHECK(std::find(right_parents_14.begin(), right_parents_14.end(), 13) ==
        right_parents_14.end());
  CHECK(std::find(right_parents_14.begin(), right_parents_14.end(), 15) !=
        right_parents_14.end());
  CHECK(std::find(right_children_14.begin(), right_children_14.end(), 11) !=
        right_children_14.end());
  CHECK_EQ(node_14.Id(), 14);
  // Check that `subsplit_to_id_` node ids were updated.
  CHECK_EQ(dag.GetDAGNodeId(node_14.GetBitset()), NodeId(14));
  // Check that `dag_edges_` node ids were updated.
  CHECK_EQ(dag.GetEdgeIdx(NodeId(15), NodeId(14)), EdgeId(11));
  // Check that `dag_edges_` edge idxs were updated.
  CHECK_EQ(dag.GetEdgeIdx(NodeId(14), NodeId(13)), EdgeId(10));
  CHECK_EQ(dag.GetEdgeIdx(NodeId(16), NodeId(13)), EdgeId(13));
  CHECK_EQ(dag.GetEdgeIdx(NodeId(11), NodeId(4)), EdgeId(25));
  // Check that `parent_to_child_range_` was updated.
  CHECK_EQ(dag.GetChildEdgeRange(node_14.GetBitset(), false).second, EdgeId(11));
  CHECK_EQ(dag.GetChildEdgeRange(dag.GetDAGNode(NodeId(16)).GetBitset(), false).first,
           EdgeId(12));
  CHECK_EQ(dag.GetChildEdgeRange(dag.GetDAGNode(NodeId(16)).GetBitset(), false).second,
           EdgeId(14));
  // Check that `topology_count_` was updated.
  CHECK_EQ(dag.TopologyCount(), prev_topology_count + 2);
}

// See diagram at https://github.com/phylovi/bito/issues/391#issuecomment-1168061363.
TEST_CASE("SubsplitDAG: Only add parent node tests") {
  const std::string fasta_path = "data/five_taxon.fasta";
  auto inst = GPInstanceOfFiles(fasta_path, "data/five_taxon_rooted_more_2.nwk");
  auto& dag = inst.GetDAG();
  // Before adding any nodes.
  size_t prev_node_count = dag.NodeCount();
  size_t prev_edge_count = dag.EdgeCountWithLeafSubsplits();
  // Add nodes 12|34 and 1|2.
  dag.AddNodePair(Bitset::Subsplit("01100", "00011"),
                  Bitset::Subsplit("01000", "00100"));
  CHECK_EQ(dag.NodeCount(), prev_node_count + 2);
  CHECK_EQ(dag.EdgeCountWithLeafSubsplits(), prev_edge_count + 5);
  // Add nodes 0|12 and 1|2 (this should just add 0|12 and associated edges).
  dag.AddNodePair(Bitset::Subsplit("10000", "01100"),
                  Bitset::Subsplit("01000", "00100"));
  // Check that the node count and edge count was updated.
  CHECK_EQ(dag.NodeCount(), prev_node_count + 3);
  CHECK_EQ(dag.EdgeCountWithLeafSubsplits(), prev_edge_count + 8);
  // Check that BuildEdgeReindexer() correctly handles left edges.
  CHECK_EQ(dag.GetChildEdgeRange(dag.GetDAGNode(NodeId(10)).GetBitset(), true).first,
           EdgeId(4));
  CHECK_EQ(dag.GetChildEdgeRange(dag.GetDAGNode(NodeId(10)).GetBitset(), true).second,
           EdgeId(6));
}

// See diagram at https://github.com/phylovi/bito/issues/391#issuecomment-1168064347.
TEST_CASE("SubsplitDAG: Only add child node tests") {
  const std::string fasta_path = "data/five_taxon.fasta";
  auto inst = GPInstanceOfFiles(fasta_path, "data/five_taxon_rooted_more_3.nwk");
  auto& dag = inst.GetDAG();
  // Before adding any nodes.
  size_t prev_node_count = dag.NodeCount();
  size_t prev_edge_count = dag.EdgeCountWithLeafSubsplits();
  // Add nodes 1|234 and 24|3 (this should just add 24|3 and associated edges).
  dag.AddNodePair(Bitset::Subsplit("01000", "00111"),
                  Bitset::Subsplit("00101", "00010"));
  // Check that the node count and edge count was updated.
  CHECK_EQ(dag.NodeCount(), prev_node_count + 1);
  CHECK_EQ(dag.EdgeCountWithLeafSubsplits(), prev_edge_count + 4);
  // Check that new child node is connected to all possible parents.
  CHECK_EQ(dag.GetChildEdgeRange(dag.GetDAGNode(NodeId(10)).GetBitset(), false).first,
           EdgeId(10));
  CHECK_EQ(dag.GetChildEdgeRange(dag.GetDAGNode(NodeId(10)).GetBitset(), false).second,
           EdgeId(12));
  CHECK_EQ(dag.GetChildEdgeRange(dag.GetDAGNode(NodeId(11)).GetBitset(), false).first,
           EdgeId(5));
  CHECK_EQ(dag.GetChildEdgeRange(dag.GetDAGNode(NodeId(11)).GetBitset(), false).second,
           EdgeId(7));
}

// Checks that parent nodes found via scan match found via map.
auto TestParentNodeIds = [](const GPDAG& dag, const Bitset& subsplit) {
  const auto [left_via_map, right_via_map] = dag.FindParentNodeIdsViaMap(subsplit);
  const auto [left_via_scan, right_via_scan] = dag.FindParentNodeIdsViaScan(subsplit);
  std::unordered_set<NodeId> left_via_map_set(left_via_map.begin(), left_via_map.end());
  std::unordered_set<NodeId> right_via_map_set(right_via_map.begin(),
                                               right_via_map.end());
  std::unordered_set<NodeId> left_via_scan_set(left_via_scan.begin(),
                                               left_via_scan.end());
  std::unordered_set<NodeId> right_via_scan_set(right_via_scan.begin(),
                                                right_via_scan.end());

  bool matches = !(left_via_map_set != left_via_scan_set or
                   right_via_map_set != right_via_scan_set);
  if (!matches) {
    std::cout << "FindParentNodeIds [FAIL_BEGIN]" << std::endl;
    std::cout << "Subsplit: " << subsplit.SubsplitToString() << std::endl;
    std::cout << "LinearSearch: " << left_via_scan_set << " " << right_via_scan_set
              << std::endl;
    std::cout << "MapSearch: " << left_via_map_set << " " << right_via_map_set
              << std::endl;

    std::cout << "via_map_set: [ ";
    for (const auto node_id : left_via_map_set) {
      std::cout << dag.GetDAGNode(node_id).GetBitset().SubsplitToString() << " ";
    }
    std::cout << "] [ ";
    for (const auto node_id : right_via_map_set) {
      std::cout << dag.GetDAGNode(node_id).GetBitset().SubsplitToString() << " ";
    }
    std::cout << "] " << std::endl;

    std::cout << "via_scan_set: [ ";
    for (const auto node_id : left_via_scan_set) {
      std::cout << dag.GetDAGNode(node_id).GetBitset().SubsplitToString() << " ";
    }
    std::cout << "] [ ";
    for (const auto node_id : right_via_scan_set) {
      std::cout << dag.GetDAGNode(node_id).GetBitset().SubsplitToString() << " ";
    }
    std::cout << "] " << std::endl;

    std::cout << "FindParentNodeIds [FAIL_END]" << std::endl;
  } else {
    // std::cout << "FindParentNodeIds [PASS]" << std::endl;
  }
  return matches;
};
// Checks that child nodes found via scan match found via map.
auto TestChildNodeIds = [](const GPDAG& dag, const Bitset& subsplit) {
  const auto [left_via_map, right_via_map] = dag.FindChildNodeIdsViaMap(subsplit);
  const auto [left_via_scan, right_via_scan] = dag.FindChildNodeIdsViaScan(subsplit);
  std::unordered_set<NodeId> left_via_map_set(left_via_map.begin(), left_via_map.end());
  std::unordered_set<NodeId> right_via_map_set(right_via_map.begin(),
                                               right_via_map.end());
  std::unordered_set<NodeId> left_via_scan_set(left_via_scan.begin(),
                                               left_via_scan.end());
  std::unordered_set<NodeId> right_via_scan_set(right_via_scan.begin(),
                                                right_via_scan.end());

  bool matches = !(left_via_map_set != left_via_scan_set or
                   right_via_map_set != right_via_scan_set);
  if (!matches) {
    std::cout << "FindChildNodeIds [FAIL_BEGIN]" << std::endl;
    std::cout << "Subsplit: " << subsplit.SubsplitToString() << std::endl;
    std::cout << "LinearSearch: " << left_via_scan_set << " " << right_via_scan_set
              << std::endl;
    std::cout << "MapSearch: " << left_via_map_set << " " << right_via_map_set
              << std::endl;

    std::cout << "via_map_set: [ ";
    for (const auto node_id : left_via_map_set) {
      std::cout << dag.GetDAGNode(node_id).GetBitset().SubsplitToString() << " ";
    }
    std::cout << "] [ ";
    for (const auto node_id : right_via_map_set) {
      std::cout << dag.GetDAGNode(node_id).GetBitset().SubsplitToString() << " ";
    }
    std::cout << "] " << std::endl;

    std::cout << "via_scan_set: [ ";
    for (const auto node_id : left_via_scan_set) {
      std::cout << dag.GetDAGNode(node_id).GetBitset().SubsplitToString() << " ";
    }
    std::cout << "] [ ";
    for (const auto node_id : right_via_scan_set) {
      std::cout << dag.GetDAGNode(node_id).GetBitset().SubsplitToString() << " ";
    }
    std::cout << "] " << std::endl;

    std::cout << "FindChildNodeIds [FAIL_END]" << std::endl;
  } else {
    // std::cout << "FindChildNodeIds [PASS]" << std::endl;
  }
  return matches;
};

// Compares adding nodes to DAG individually vs adding multiple nodes. Additionally,
// tests that adjacent nodes from acquired via map lookup match those acquired via
// linear scan.
TEST_CASE("SubsplitDAG: Add Multiple Nodes") {
  const std::string fasta_path = "data/six_taxon.fasta";
  const std::string newick_path = "data/six_taxon_rooted_simple.nwk";
  // Instance that will be unaltered.
  auto inst = GPInstanceOfFiles(fasta_path, newick_path);
  inst.MakeNNIEngine();
  NNIEngine& nni_engine = inst.GetNNIEngine();
  GPDAG dag1 = inst.GetDAG();
  GPDAG dag2 = inst.GetDAG();
  nni_engine.SyncAdjacentNNIsWithDAG();

  // Check unaltered DAG nodes match via map and via linear scan.
  for (const auto node_id : dag1.LeafwardNodeTraversalTrace(true)) {
    const auto subsplit = dag1.GetDAGNodeBitset(node_id);
    CHECK_MESSAGE(TestChildNodeIds(dag1, subsplit),
                  "Child nodes found by map lookup do not match those found by linear "
                  "scan (before adding nodes).");
    CHECK_MESSAGE(TestParentNodeIds(dag1, subsplit),
                  "Parent nodes found by map lookup do not match those found by linear "
                  "scan (before adding nodes).");
  }

  // Check DAG nodes match after adding nodes individually.
  for (const auto nni : nni_engine.GetAdjacentNNIs()) {
    auto mods1 = dag1.AddNodePair(nni);
    for (const auto node_id : dag1.LeafwardNodeTraversalTrace(true)) {
      const auto subsplit = dag1.GetDAGNodeBitset(node_id);
      CHECK_MESSAGE(TestChildNodeIds(dag1, subsplit),
                    "Child nodes found by map lookup do not match those found by "
                    "linear scan (after adding nodes).");
      CHECK_MESSAGE(TestParentNodeIds(dag1, subsplit),
                    "Parent nodes found by map lookup do not match those found by "
                    "linear scan (after adding nodes).");
    }
  }

  CHECK_MESSAGE(
      dag1 != dag2,
      "DAG with nodes added individually incorrectly matches the DAG with nodes "
      "added collectively (before adding nodes).");
  // Check DAGs match by adding nodes individually vs all-at-once.
  std::vector<std::pair<Bitset, Bitset>> node_subsplit_pairs;
  for (const auto nni : nni_engine.GetAdjacentNNIs()) {
    node_subsplit_pairs.push_back({nni.GetParent(), nni.GetChild()});
  }
  dag2.AddNodes(node_subsplit_pairs);
  CHECK_MESSAGE(dag1 == dag2,
                "DAG with nodes added individually does not match the DAG with nodes "
                "added collectively.");
}

// Compares adding nodes to DAG individually vs adding multiple nodes. Additionally,
// tests that adjacent nodes from acquired via map lookup match those acquired via
// linear scan.
TEST_CASE("SubsplitDAG: Graft Multiple Nodes") {
  const std::string fasta_path = "data/six_taxon.fasta";
  const std::string newick_path = "data/six_taxon_rooted_simple.nwk";
  // Instance that will be unaltered.
  auto inst = GPInstanceOfFiles(fasta_path, newick_path);
  inst.MakeNNIEngine();
  NNIEngine& nni_engine = inst.GetNNIEngine();
  GPDAG dag1 = inst.GetDAG();
  GraftDAG& graft_dag = nni_engine.GetGraftDAG();
  nni_engine.SyncAdjacentNNIsWithDAG();

  // Check DAG nodes match after grafting nodes individually.
  for (const auto nni : nni_engine.GetAdjacentNNIs()) {
    auto mods1 = graft_dag.AddNodePair(nni);
    for (const auto node_id : dag1.LeafwardNodeTraversalTrace(true)) {
      const auto subsplit = dag1.GetDAGNodeBitset(node_id);
      CHECK_MESSAGE(TestChildNodeIds(dag1, subsplit),
                    "Child nodes found by map lookup do not match those found by "
                    "linear scan (after adding graft).");
      CHECK_MESSAGE(TestParentNodeIds(dag1, subsplit),
                    "Parent nodes found by map lookup do not match those found by "
                    "linear scan (after adding graft).");
    }
    graft_dag.RemoveAllGrafts();
    for (const auto node_id : dag1.LeafwardNodeTraversalTrace(true)) {
      const auto subsplit = dag1.GetDAGNodeBitset(node_id);
      CHECK_MESSAGE(TestChildNodeIds(dag1, subsplit),
                    "Child nodes found by map lookup do not match those found by "
                    "linear scan (after removing graft).");
      CHECK_MESSAGE(TestParentNodeIds(dag1, subsplit),
                    "Parent nodes found by map lookup do not match those found by "
                    "linear scan (after removing graft).");
    }
  }
}

// Tests DAG after modifying
TEST_CASE("SubsplitDAG: Add Multiple Edges") {
  const std::string fasta_path = "data/six_taxon.fasta";
  const std::string newick_path = "data/six_taxon_rooted_simple.nwk";
  // Instance that will be unaltered.
  auto inst = GPInstanceOfFiles(fasta_path, newick_path);
  inst.MakeNNIEngine();
  NNIEngine& nni_engine = inst.GetNNIEngine();
  GPDAG dag2 = inst.GetDAG();
  nni_engine.SyncAdjacentNNIsWithDAG();

  for (const auto nni : nni_engine.GetAdjacentNNIs()) {
    GPDAG dag1 = inst.GetDAG();
    // Get connecting nodes.
    const auto grandparent_nodeid = dag1.FindFirstParentNodeId(nni.GetParent());
    const auto sister_nodeid =
        dag1.FindFirstChildNodeId(nni.GetParent(), nni.WhichCladeIsSister());
    const auto leftchild_nodeid =
        dag1.FindFirstChildNodeId(nni.GetChild(), SubsplitClade::Left);
    const auto rightchild_nodeid =
        dag1.FindFirstChildNodeId(nni.GetChild(), SubsplitClade::Right);
    std::vector<NodeId> node_ids{
        {grandparent_nodeid, sister_nodeid, leftchild_nodeid, rightchild_nodeid}};
    // Get PCSP Bitsets.
    std::vector<Bitset> pcsps;
    pcsps.push_back(
        Bitset::PCSP(dag1.GetDAGNodeBitset(grandparent_nodeid), nni.GetParent()));
    pcsps.push_back(
        Bitset::PCSP(nni.GetParent(), dag1.GetDAGNodeBitset(sister_nodeid)));
    pcsps.push_back(nni.GetCentralEdgePCSP());
    pcsps.push_back(
        Bitset::PCSP(nni.GetChild(), dag1.GetDAGNodeBitset(leftchild_nodeid)));
    pcsps.push_back(
        Bitset::PCSP(nni.GetChild(), dag1.GetDAGNodeBitset(rightchild_nodeid)));
    dag1.AddEdges(pcsps);
    for (const auto pcsp : pcsps) {
      CHECK_MESSAGE(dag1.ContainsEdge(pcsp), "DAG does not contain added edge.");
    }
    CHECK_MESSAGE(
        dag1.EdgeCountWithLeafSubsplits() == dag2.EdgeCountWithLeafSubsplits() + 5,
        "DAG does not contain proper number of edges after adding edges.");
  }
}

// ** NNIEngine tests **

using NodeMap = std::unordered_map<NodeId, NodeId>;
using EdgeMap = std::unordered_map<EdgeId, EdgeId>;
using TreeIdMap = std::unordered_map<TreeId, RootedTree>;
using TreeEdgeMap = std::unordered_map<EdgeId, TreeId>;
using EdgeScoreMap = std::unordered_map<EdgeId, double>;
using TreeScoreMap = std::unordered_map<TreeId, double>;
using NNIScoreMap = std::map<NNIOperation, double>;
using BranchLengths = EigenVectorXd;
using NNIBranchLengthsMap = std::map<NNIOperation, EigenVectorXd>;
using BranchMap = DAGBranchHandler::BranchLengthMap;
using NNIBranchMapMap = std::map<NNIOperation, DAGBranchHandler::BranchLengthMap>;

// Builds a mapping from node and edge elements from pre-DAG to post-DAG.  If pre-DAG
// element does not exist in post-DAG, return NoId.
std::pair<NodeMap, EdgeMap> BuildNodeAndEdgeMapsFromPreDAGToPostDAG(GPDAG& pre_dag,
                                                                    GPDAG& post_dag) {
  NodeMap node_map;
  for (const auto& bitset : pre_dag.BuildSetOfNodeBitsets()) {
    if (bitset.SubsplitIsUCA()) {
      continue;
    }
    const auto pre_id = pre_dag.GetDAGNodeId(bitset);
    auto post_id = NodeId(NodeId::NoId);
    if (post_dag.ContainsNode(bitset)) {
      post_id = post_dag.GetDAGNodeId(bitset);
    }
    node_map[pre_id] = post_id;
  }

  EdgeMap edge_map;
  for (const auto& bitset : pre_dag.BuildSetOfEdgeBitsets()) {
    const auto pre_id = pre_dag.GetEdgeIdx(bitset);
    auto post_id = EdgeId(EdgeId::NoId);
    if (post_dag.ContainsEdge(bitset)) {
      post_id = post_dag.GetEdgeIdx(bitset);
    }
    edge_map[pre_id] = post_id;
  }

  return std::make_pair(node_map, edge_map);
}

// Tests that reindexers match the remapped node_ids and edge_idxs after AddNodePair.
TEST_CASE("NNIEngine: Reindexers for AddNodePair") {
  const std::string fasta_path = "data/five_taxon.fasta";
  const std::string newick_path = "data/five_taxon_rooted_more_2.nwk";
  auto pre_inst = GPInstanceOfFiles(fasta_path, newick_path);
  auto& pre_dag = pre_inst.GetDAG();
  auto inst = GPInstanceOfFiles(fasta_path, newick_path);
  auto& dag = inst.GetDAG();
  inst.MakeNNIEngine();
  auto& nni_engine = inst.GetNNIEngine();
  nni_engine.SyncAdjacentNNIsWithDAG();

  for (const auto& nni : nni_engine.GetAdjacentNNIs()) {
    auto mods = dag.AddNodePair(nni);
    for (NodeId old_idx = NodeId(0); old_idx < pre_dag.NodeCount(); old_idx++) {
      NodeId new_idx =
          NodeId(mods.node_reindexer.GetNewIndexByOldIndex(old_idx.value_));
      Bitset old_node = pre_dag.GetDAGNode(old_idx).GetBitset();
      Bitset new_node = dag.GetDAGNode(new_idx).GetBitset();
      CHECK_EQ(old_node, new_node);
    }
    for (EdgeId old_idx = EdgeId(0); old_idx < pre_dag.EdgeCount(); old_idx++) {
      EdgeId new_idx =
          EdgeId(mods.edge_reindexer.GetNewIndexByOldIndex(old_idx.value_));
      Bitset old_parent =
          pre_dag.GetDAGNode(NodeId(pre_dag.GetDAGEdge(old_idx).GetParent()))
              .GetBitset();
      Bitset old_child =
          pre_dag.GetDAGNode(NodeId(pre_dag.GetDAGEdge(old_idx).GetChild()))
              .GetBitset();
      Bitset new_parent =
          dag.GetDAGNode(NodeId(dag.GetDAGEdge(new_idx).GetParent())).GetBitset();
      Bitset new_child =
          dag.GetDAGNode(NodeId(dag.GetDAGEdge(new_idx).GetChild())).GetBitset();
      CHECK_EQ(old_parent, new_parent);
      CHECK_EQ(old_child, new_child);
    }
    pre_dag.AddNodePair(nni);
  }
}

// This test builds a DAG, tests if engine generates the same set of adjacent NNIs and
// manually created set. Then adds a node pair to DAG, and tests if engine updates
// correctly.
TEST_CASE("NNIEngine: Adjacent NNI Maintenance") {
  // Simple DAG that contains a shared edge, internal leafward fork, and an internal
  // rootward fork.
  const std::string fasta_path = "data/six_taxon.fasta";
  auto inst = GPInstanceOfFiles(fasta_path, "data/six_taxon_rooted_simple.nwk");
  GPDAG& dag = inst.GetDAG();
  GPEngine& gp_engine = inst.GetGPEngine();

  NNISet correct_adjacent_nnis;

  auto nni_engine = NNIEngine(dag, &gp_engine);
  auto nni_engine_2 = NNIEngine(dag, &gp_engine);

  // Build adjacent NNIs from current DAG state.
  nni_engine.SyncAdjacentNNIsWithDAG();

  // Functions for quick manual insertion/removal for Correct NNI Set.
  auto InsertNNI = [&correct_adjacent_nnis](Bitset parent, Bitset child) {
    auto nni = NNIOperation(parent, child);
    correct_adjacent_nnis.insert(nni);
  };
  auto RemoveNNI = [&correct_adjacent_nnis](Bitset parent, Bitset child) {
    auto nni = NNIOperation(parent, child);
    correct_adjacent_nnis.erase(nni);
  };

  // For images and notes describing this part of the test case, see
  // https://github.com/phylovi/bito/pull/366#issuecomment-920454401
  // Add NNIs for edge 4 to NNI Set.
  InsertNNI(Bitset::Subsplit("010000", "101111"),   //  (parent)-(child)
            Bitset::Subsplit("100000", "001111"));  // (1|02345)-(0|2345)
  InsertNNI(Bitset::Subsplit("100000", "011111"),
            Bitset::Subsplit("010000", "001111"));  // (0|12345)-(1|2345)
  // Add NNIs for edge 6 to NNI Set.
  InsertNNI(Bitset::Subsplit("001000", "110111"),
            Bitset::Subsplit("110000", "000111"));  // (2|01345)-(01|345)
  InsertNNI(Bitset::Subsplit("000111", "111000"),
            Bitset::Subsplit("110000", "001000"));  // (345|012)-(01|2)
  // Add NNIs for edge 7 to NNI Set.
  InsertNNI(Bitset::Subsplit("000001", "111110"),
            Bitset::Subsplit("110000", "001110"));  // (5|01234)-(01|234)
  InsertNNI(Bitset::Subsplit("001110", "110001"),
            Bitset::Subsplit("110000", "000001"));  // (234|015)-(01|5)
  // Add NNIs for edge 2 to NNI Set.
  InsertNNI(Bitset::Subsplit("000110", "001001"),
            Bitset::Subsplit("001000", "000001"));  // (34|25)-(2|5)
  // No NNIs to add for edge 5 to NNI Set (see notes).
  // Add NNIs for edge 3 to NNI Set.
  InsertNNI(Bitset::Subsplit("000100", "001010"),
            Bitset::Subsplit("001000", "000010"));  // (3|24)-(2|4)
  InsertNNI(Bitset::Subsplit("000010", "001100"),
            Bitset::Subsplit("001000", "000100"));  // (4|23)-(2|3)
  // Add NNIs for edge 1 to NNI Set.
  InsertNNI(Bitset::Subsplit("000010", "000101"),
            Bitset::Subsplit("000100", "000001"));  // (4|35)-(3|5)
  InsertNNI(Bitset::Subsplit("000100", "000011"),
            Bitset::Subsplit("000010", "000001"));  // (3|45)-(4|5)

  // Check that `BuildNNISet()` added correct set of nnis.
  auto adjacent_nnis = nni_engine.GetAdjacentNNIs();
  CHECK_EQ(adjacent_nnis, correct_adjacent_nnis);

  // Now we add a node pair to DAG so we can check UpdateAdjacentNNIsAfterAddNodePair.
  // see https://github.com/phylovi/bito/pull/366#issuecomment-922781415
  NNIOperation nni_to_add(Bitset::Subsplit("000110", "001001"),
                          Bitset::Subsplit("001000", "000001"));  // (34|25)-(2|5)
  dag.AddNodePair(nni_to_add);

  // Update NNI.
  nni_engine.UpdateAdjacentNNIsAfterDAGAddNodePair(nni_to_add);
  // Add parents of parent (edge 8) to NNI Set.
  InsertNNI(Bitset::Subsplit("001001", "110110"),
            Bitset::Subsplit("110000", "000110"));  // (25|0134)-(01|34)
  InsertNNI(Bitset::Subsplit("000110", "111001"),
            Bitset::Subsplit("110000", "001001"));  // (34|0125)-(01|25)
  // Add children of parent (edge 19) to NNI Set.
  InsertNNI(Bitset::Subsplit("000100", "001011"),
            Bitset::Subsplit("001001", "000010"));  // (3|245)-(25|4)
  InsertNNI(Bitset::Subsplit("000010", "001101"),
            Bitset::Subsplit("001001", "000100"));  // (4|235)-(25|3)
  // No parents of child (edge 20) to add to NNI Set (see notes).
  // These should not be equal, as it has not yet removed the pair just added to DAG.
  CHECK_NE(nni_engine.GetAdjacentNNIs(), correct_adjacent_nnis);
  // Remove NNI added to DAG from NNI Set.
  RemoveNNI(nni_to_add.parent_, nni_to_add.child_);
  // Check that `UpdateAdjacentNNIsAfterAddNodePair()` updated correctly.
  CHECK_EQ(nni_engine.GetAdjacentNNIs(), correct_adjacent_nnis);

  // Build NNI Set from current DAG state from scratch.
  nni_engine_2.SyncAdjacentNNIsWithDAG();
  CHECK_EQ(nni_engine_2.GetAdjacentNNIs(), correct_adjacent_nnis);
}

// Tests DAG equality after adding different NNIs and built from different taxon
// orderings. Test described at:
// https://github.com/phylovi/bito/pull/377#issuecomment-1035410447
TEST_CASE("NNIEngine: Add NNI Test") {
  // Fasta contains simple sequences for four taxa: x0,x1,x2,x3.
  const std::string fasta_path = "data/four_taxon.fasta";
  // dag_A_1 is a DAG that contains pair_1.
  auto inst_A_1 =
      GPInstanceOfFiles(fasta_path, "data/four_taxon_simple_before_nni_1.nwk",
                        "_ignore/mmapped_pv_A_1.data");
  GPDAG& dag_A_1 = inst_A_1.GetDAG();
  // dag_A_2 is a DAG that contains pair_2.
  auto inst_A_2 =
      GPInstanceOfFiles(fasta_path, "data/four_taxon_simple_before_nni_2.nwk",
                        "_ignore/mmapped_pv_A_2.data");
  GPDAG& dag_A_2 = inst_A_2.GetDAG();
  // dag_A_2b is a DAG that contains pair_2 with a different taxon mapping.
  auto inst_A_2b =
      GPInstanceOfFiles(fasta_path, "data/four_taxon_simple_before_nni_2b.nwk",
                        "_ignore/mmapped_pv_A_2.data");
  GPDAG& dag_A_2b = inst_A_2b.GetDAG();
  // dag_B is a DAG containing dag_A_1 after adding node pair_2.
  auto inst_B = GPInstanceOfFiles(fasta_path, "data/four_taxon_simple_after_nni.nwk",
                                  "_ignore/mmapped_pv_B.data");
  GPDAG& dag_B = inst_B.GetDAG();
  // pair_1: NNI pair missing from dag_A_1.
  NNIOperation pair_1(Bitset::Subsplit("0110", "0001"),   // 12|3
                      Bitset::Subsplit("0100", "0010"));  //  1|2
  // pair_2: NNI pair missing from dag_A_2.
  NNIOperation pair_2(Bitset::Subsplit("0110", "0001"),   // 12|3
                      Bitset::Subsplit("0100", "0010"));  //  1|2
  // pair_2b: NNI pair missing from dag_A_2b.
  NNIOperation pair_2b(Bitset::Subsplit("0100", "0011"),   //  1|23
                       Bitset::Subsplit("0010", "0001"));  //  2|3
  // Before adding missing NNIs, dag_A_2 variants are equal, but dag_A_1 and dag_A_2 are
  // different.
  CHECK_EQ(dag_A_1, dag_A_1);
  CHECK_EQ(dag_A_2, dag_A_2b);
  CHECK_NE(dag_A_1, dag_A_2);
  // Add missing NNIs.
  dag_A_1.AddNodePair(pair_1);
  dag_A_2.AddNodePair(pair_2);
  dag_A_2b.AddNodePair(pair_2b);
  // After adding missing NNIs, all DAGs are equal to dag_B.
  CHECK_EQ(dag_A_1, dag_B);
  CHECK_EQ(dag_A_2, dag_B);
  CHECK_EQ(dag_A_2b, dag_B);
}

// Starts with a DAG built from a single tree. Iteratively finds all adjacent NNIs and
// adds them to the DAG, until there are no more adjacent NNIs to DAG.
// (1) Tests that resulting DAG is equal to the complete DAG, containing all possible
// subsplits. (2) Reruns with "include rootsplit" option off.  Checks that
// resulting DAG contains only edges reachable from the initial rootsplit.
TEST_CASE("NNIEngine: Build Complete DAG by Adding NNIs (include/exclude rootsplit)") {
  auto BuildCompleteDAGSubsplits = [](const size_t taxon_count) {
    auto BuildSubsplitsFromAllTaxonAssignment =
        [](std::set<Bitset>& all_subsplits, Bitset& subsplit, const size_t i,
           auto&& BuildSubsplitsFromAllTaxonAssignment) {
          if (i == subsplit.SubsplitGetCladeSize()) {
            auto subsplit_out =
                Bitset::Subsplit(subsplit.SubsplitGetClade(SubsplitClade::Left),
                                 subsplit.SubsplitGetClade(SubsplitClade::Right));
            if (subsplit_out.SubsplitGetClade(SubsplitClade::Right).None() ||
                subsplit_out.SubsplitGetClade(SubsplitClade::Left).None()) {
              return;
            }
            all_subsplits.insert(subsplit_out);
            return;
          }
          for (size_t j = 0; j < 3; j++) {
            if (j == 0) {
              subsplit.set(i, true);
              subsplit.set(subsplit.SubsplitGetCladeSize() + i, false);
            } else if (j == 1) {
              subsplit.set(i, false);
              subsplit.set(subsplit.SubsplitGetCladeSize() + i, true);
            } else {
              subsplit.set(i, false);
              subsplit.set(subsplit.SubsplitGetCladeSize() + i, false);
            }
            BuildSubsplitsFromAllTaxonAssignment(all_subsplits, subsplit, i + 1,
                                                 BuildSubsplitsFromAllTaxonAssignment);
          }
        };

    std::set<Bitset> subsplits;
    Bitset subsplit(taxon_count * 2, false);
    BuildSubsplitsFromAllTaxonAssignment(subsplits, subsplit, 0,
                                         BuildSubsplitsFromAllTaxonAssignment);
    for (size_t i = 0; i < taxon_count; i++) {
      subsplits.insert(
          Bitset::LeafSubsplitOfNonemptyClade(Bitset::Singleton(taxon_count, i)));
    }
    subsplits.insert(Bitset::UCASubsplitOfTaxonCount(taxon_count));

    return subsplits;
  };

  auto TestCompleteDAG = [&BuildCompleteDAGSubsplits](
                             const std::string& fasta_path,
                             const std::string& newick_path,
                             const bool include_rootsplits = true) {
    auto inst = GPInstanceOfFiles(fasta_path, newick_path);
    auto& dag = inst.GetDAG();
    std::set<Bitset> rootsplit_subsplits;
    for (const auto& node_id : dag.GetRootsplitNodeIds()) {
      const auto& subsplit = dag.GetDAGNodeBitset(node_id);
      rootsplit_subsplits.insert(subsplit);
    }

    inst.MakeNNIEngine();
    auto& nniengine = inst.GetNNIEngine();
    nniengine.SetIncludeRootsplitNNIs(include_rootsplits);
    nniengine.SyncAdjacentNNIsWithDAG();
    while (nniengine.GetAdjacentNNICount() > 0) {
      for (const auto& nni : nniengine.GetAdjacentNNIs()) {
        dag.AddNodePair(nni.GetParent(), nni.GetChild());
      }
      nniengine.SyncAdjacentNNIsWithDAG();
    }

    auto subsplits = BuildCompleteDAGSubsplits(dag.TaxonCount());
    if (!include_rootsplits) {
      std::set<Bitset> tmp_rootsplit_subsplits;
      for (const auto& node_id : dag.GetRootsplitNodeIds()) {
        const auto& subsplit = dag.GetDAGNodeBitset(node_id);
        tmp_rootsplit_subsplits.insert(subsplit);
      }
      CHECK_MESSAGE(
          rootsplit_subsplits == tmp_rootsplit_subsplits,
          "Rootsplit NNIs were added when NNIEngine flagged to exclude rootsplits.");

      std::set<Bitset> tmp_subsplits;
      for (const auto& ancestor : rootsplit_subsplits) {
        tmp_subsplits.insert(ancestor);
      }
      for (const auto& descendant : subsplits) {
        for (const auto& ancestor : rootsplit_subsplits) {
          if (Bitset::SubsplitIsAncestorDescendantPair(ancestor, descendant,
                                                       SubsplitClade::Left) ||
              Bitset::SubsplitIsAncestorDescendantPair(ancestor, descendant,
                                                       SubsplitClade::Right)) {
            tmp_subsplits.insert(descendant);
            break;
          }
        }
      }
      tmp_subsplits.insert(Bitset::UCASubsplitOfTaxonCount(dag.TaxonCount()));
      subsplits = tmp_subsplits;
    }

    bool contains_all_subsplits = true;
    for (const auto& subsplit : subsplits) {
      contains_all_subsplits &= dag.ContainsNode(subsplit);
      if (!contains_all_subsplits) {
        std::cout << "Missing subsplit from complete DAG: "
                  << subsplit.SubsplitToString() << std::endl;
        break;
      }
    }
    contains_all_subsplits &= (subsplits.size() == dag.NodeCount());
    CHECK_MESSAGE(subsplits.size() == dag.NodeCount(),
                  "DAG node count is not equal to number of nodes in Complete DAG.");
    return contains_all_subsplits;
  };

  const std::string fasta_path_0 = "data/hello.fasta";
  const std::string newick_path_0 = "data/hello_rooted_diff_branches.nwk";
  CHECK_MESSAGE(TestCompleteDAG(fasta_path_0, newick_path_0, true),
                "Complete DAG test for Hello.");
  CHECK_MESSAGE(TestCompleteDAG(fasta_path_0, newick_path_0, false),
                "Complete DAG test for Hello when excluding rootsplits.");

  const std::string fasta_path_1 = "data/five_taxon.fasta";
  const std::string newick_path_1 = "data/five_taxon_trees_3_4_diff_branches.nwk";
  CHECK_MESSAGE(TestCompleteDAG(fasta_path_1, newick_path_1, true),
                "Complete DAG test for Five Taxon.");
  CHECK_MESSAGE(TestCompleteDAG(fasta_path_1, newick_path_1, false),
                "Complete DAG test for Five Taxon when excluding rootsplits.");
}

// Access ith NNI from NNI set.
NNIOperation GetWhichNNIFromSet(const NNISet& nni_set, const size_t which_nni) {
  auto nni_set_ptr = nni_set.begin();
  for (size_t i = 0; i < which_nni; i++) {
    nni_set_ptr++;
  }
  return *nni_set_ptr;
};

// This compares DAGs after adding NNIs to SubsplitDAG vs GraftDAG.
// Also tests that Adding all node pairs to DAG gives proper result.
TEST_CASE("NNIEngine: GraftDAG") {
  // Simple DAG that contains a shared edge, internal leafward fork, and an internal
  // rootward fork.
  const std::string fasta_path = "data/six_taxon.fasta";
  const std::string newick_path = "data/six_taxon_rooted_simple.nwk";
  // Instance that will be unaltered.
  auto pre_inst = GPInstanceOfFiles(fasta_path, newick_path);
  GPDAG& pre_dag = pre_inst.GetDAG();
  // Instance that is used by grafted DAG.
  auto graft_inst = GPInstanceOfFiles(fasta_path, newick_path);
  graft_inst.MakeNNIEngine();
  NNIEngine& nni_engine = graft_inst.GetNNIEngine();
  GPDAG& host_dag = graft_inst.GetDAG();
  GraftDAG& graft_dag = nni_engine.GetGraftDAG();
  // Find NNIs of DAG.
  nni_engine.SyncAdjacentNNIsWithDAG();
  size_t nni_count = nni_engine.GetAdjacentNNICount();

  // TEST #0:
  // Check DAG and GraftDAG equal before adding any NNIs.
  CHECK_MESSAGE(GraftDAG::CompareToDAG(graft_dag, pre_dag) == 0,
                "GraftDAG not equal to DAG before AddNodePair.");
  // TEST #1:
  // Add each NNI pair to GraftDAG, then check that expected nodes and edges are present
  // in DAG.
  for (size_t i = 0; i < nni_count; i++) {
    auto nni = GetWhichNNIFromSet(nni_engine.GetAdjacentNNIs(), i);
    graft_dag.AddNodePair(nni.parent_, nni.child_);
    CHECK_MESSAGE(graft_dag.ContainsNode(nni.parent_),
                  "Cannot find parent node in GraftDAG.");
    CHECK_MESSAGE(graft_dag.ContainsNode(nni.child_),
                  "Cannot find child node in GraftDAG.");
    if (!(graft_dag.ContainsNode(nni.parent_) || graft_dag.ContainsNode(nni.child_))) {
      graft_dag.RemoveAllGrafts();
      continue;
    }
    size_t neighbor_count = 0;
    for (const auto is_parent : {true, false}) {
      const auto& subsplit = (is_parent ? nni.parent_ : nni.child_);
      const auto node_id = graft_dag.GetDAGNodeId(subsplit);
      const auto& node = graft_dag.GetDAGNode(node_id);
      for (const auto direction : {true, false}) {
        for (const auto clade : {true, false}) {
          const auto neighbors = node.GetEdge(direction, clade);
          for (const auto neighbor_id : neighbors) {
            const auto parent_id = NodeId(is_parent ? node_id : neighbor_id);
            const auto child_id = NodeId(is_parent ? neighbor_id : node_id);
            CHECK_MESSAGE(graft_dag.ContainsEdge(parent_id, child_id),
                          "Cannot find edge in GraftDAG.");
            CHECK_MESSAGE(graft_dag.ContainsEdge(child_id, parent_id),
                          "Cannot find edge in GraftDAG.");
            neighbor_count++;
          }
        }
      }
    }
    if (neighbor_count < 6) {
      std::cout << "Too few neighbors to NNI: " << i << " " << nni << std::endl;
    }
    CHECK_MESSAGE(neighbor_count >= 6, "There were too few neighbors to NNI.");
  }
  graft_dag.RemoveAllGrafts();
  // TEST #2:
  // Test DAGs equivalent after adding NNI individually.
  for (size_t i = 0; i < nni_count; i++) {
    auto nni = GetWhichNNIFromSet(nni_engine.GetAdjacentNNIs(), i);
    // New temp DAG.
    auto inst = GPInstanceOfFiles(fasta_path, newick_path);
    GPDAG& dag = inst.GetDAG();
    // Add NNI to DAG and GraftDAG and compare results.
    dag.AddNodePair(nni.parent_, nni.child_);
    graft_dag.AddNodePair(nni.parent_, nni.child_);
    CHECK_MESSAGE(GraftDAG::CompareToDAG(graft_dag, dag) == 0,
                  "GraftDAG not equal to DAG after adding NNIs.");
    // Clear NNIs from GraftDAG and compare to initial DAG.
    graft_dag.RemoveAllGrafts();
    CHECK_MESSAGE(GraftDAG::CompareToDAG(graft_dag, pre_dag) == 0,
                  "GraftDAG not equal to initial DAG after removing all NNIs.");
  }
  // TEST #3:
  // Test DAGs not equivalent when adding different NNIs.
  {
    auto nni_1 = GetWhichNNIFromSet(nni_engine.GetAdjacentNNIs(), 0);
    auto nni_2 = GetWhichNNIFromSet(nni_engine.GetAdjacentNNIs(), 1);
    // New temp DAG.
    auto inst = GPInstanceOfFiles(fasta_path, newick_path);
    GPDAG& dag = inst.GetDAG();
    // Add NNI to DAG and GraftDAG and compare results.
    dag.AddNodePair(nni_1.parent_, nni_1.child_);
    graft_dag.AddNodePair(nni_2.parent_, nni_2.child_);
    CHECK_MESSAGE(GraftDAG::CompareToDAG(graft_dag, dag) != 0,
                  "GraftDAG is equal to DAG after adding different NNIs.");
  }

  // TEST #4:
  // Modify GraftDAG, clear GraftDAG, modify DAG, then modify GraftDAG again.
  for (size_t i = 0; i < nni_count; i++) {
    auto nni = GetWhichNNIFromSet(nni_engine.GetAdjacentNNIs(), i);
    graft_dag.AddNodePair(nni.parent_, nni.child_);

    CHECK_MESSAGE(
        (graft_dag.ContainsNode(nni.parent_) && graft_dag.ContainsNode(nni.child_)),
        "Graft DAG does not contain added nodes by the Graft DAG.");
    graft_dag.RemoveAllGrafts();
    host_dag.AddNodePair(nni.parent_, nni.child_);
    CHECK_MESSAGE(
        (host_dag.ContainsNode(nni.parent_) && host_dag.ContainsNode(nni.child_)),
        "Host DAG does not contain added nodes added by the Host DAG.");
    CHECK_MESSAGE(
        (graft_dag.ContainsNode(nni.parent_) && graft_dag.ContainsNode(nni.child_)),
        "Graft DAG does not contain added nodes added by the Graft DAG.");
    CHECK_MESSAGE(
        graft_dag.GetDAGNodeId(nni.parent_) == host_dag.GetDAGNodeId(nni.parent_),
        "Graft DAG NodeId does not match Host DAG NodeId.");
  }
}

// Initialize GPInstance, make GPEngine, DAG, and NNIEngine.
// Perform initial run of GP optimization.
void GPInstanceRunGP(GPInstance& inst, const bool do_optimize_branch_lengths = true,
                     const bool do_fixed_branch_lengths = false,
                     const bool do_reinit_priors = true) {
  if (do_fixed_branch_lengths) {
    inst.GetGPEngine().SetBranchLengthsToDefault();
  }
  if (do_optimize_branch_lengths) {
    inst.EstimateBranchLengths(0.0001, 100, true);
  }
  if (do_reinit_priors) {
    inst.ReinitializePriors();
  }
  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  inst.ComputeMarginalLikelihood();
}

// Adds NNIs to DAG, then resizes and reindexes GPEngine, then checks that the same
// node and edge bitsets correspond to the same PLVs and branch lengths before and after
// AddNodePair.
TEST_CASE("NNIEngine: Resize and Reindex GPEngine after AddNodePair") {
  // Check that two GPInstances produce the same results after GP run.
  auto CheckGPEngineRun = [](GPInstance& inst, GPInstance& pre_inst) {
    bool passes_gp_run = true;
    inst.EstimateBranchLengths(0.0001, 100, true);
    inst.PopulatePLVs();
    inst.ComputeLikelihoods();
    inst.ComputeMarginalLikelihood();
    auto likelihoods = inst.GetGPEngine().GetPerGPCSPLogLikelihoods();
    pre_inst.EstimateBranchLengths(0.0001, 100, true);
    pre_inst.PopulatePLVs();
    pre_inst.ComputeLikelihoods();
    pre_inst.ComputeMarginalLikelihood();
    auto pre_likelihoods = inst.GetGPEngine().GetPerGPCSPLogLikelihoods();
    if (!VectorXdEquality(likelihoods, pre_likelihoods, 1e-3)) {
      return false;
    };
    return passes_gp_run;
  };
  // Check that GPEngine resized and reindexed after DAG modifications.
  auto CheckGPEngineResizeAndReindex = [](GPDAG& dag, GPEngine& gpengine,
                                          GPDAG& pre_dag, GPEngine& pre_gpengine) {
    bool passes_resized = true;
    bool passes_plv_reindexed = true;
    bool passes_branch_reindexed = true;
    // Check resizing GPEngine properly.
    passes_resized &= (gpengine.GetNodeCount() == dag.NodeCountWithoutDAGRoot());
    passes_resized &= (gpengine.GetGPCSPCount() == dag.EdgeCountWithLeafSubsplits());
    // Check that elements reindexed properly.
    const auto& [node_map, edge_map] =
        BuildNodeAndEdgeMapsFromPreDAGToPostDAG(pre_dag, dag);
    // Check PLVs reindexed properly.
    for (const auto& [pre_node_id, node_id] : node_map) {
      for (const auto plv_type : PLVTypeEnum::Iterator()) {
        const auto& plv_a = gpengine.GetPLVHandler().GetPV(plv_type, node_id);
        const auto& plv_b = pre_gpengine.GetPLVHandler().GetPV(plv_type, pre_node_id);
        auto max_diff = PLVNodeHandler::MaxDifference(plv_a, plv_b);
        if (max_diff > 1e-3) {
          passes_plv_reindexed = false;
        }
      }
    }
    // Check branch length reindexed properly.
    auto& pre_branch_lengths = pre_gpengine.GetBranchLengthHandler();
    auto& branch_lengths = gpengine.GetBranchLengthHandler();
    for (const auto& [pre_edge_id, edge_id] : edge_map) {
      const auto branch_a = branch_lengths.Get(edge_id);
      const auto branch_b = pre_branch_lengths.Get(pre_edge_id);
      if (abs(branch_a - branch_b) > 1e-3) {
        passes_branch_reindexed = false;
      }
    }
    return passes_resized && passes_plv_reindexed && passes_branch_reindexed;
  };
  // Test that adds nodes to DAG, resizes and reindexes GPEngine and checks that
  // GPEngine reindexed correctly.
  auto ResizeAndReindexGPEngineTest = [&CheckGPEngineResizeAndReindex,
                                       &CheckGPEngineRun](
                                          const size_t nni_add_limit,
                                          const size_t test_after_every,
                                          const bool skip_reindexing,
                                          const bool perform_resize_unmodded_test) {
    BoolVector test_array;
    bool test_passes = true;
    const std::string fasta_path = "data/hotstart.fasta";
    const std::string newick_path = "data/hotstart_bootstrap_sample.nwk";
    // Instance that will not be modified.
    auto pre_inst =
        GPInstanceOfFiles(fasta_path, newick_path, "_ignore/mmap_plv_A.data");
    GPInstanceRunGP(pre_inst);
    GPDAG& pre_dag = pre_inst.GetDAG();
    GPEngine& pre_gpengine = pre_inst.GetGPEngine();
    // Instance that will have DAG and GPEngine modified.
    auto inst = GPInstanceOfFiles(fasta_path, newick_path, "_ignore/mmap_plv_C.data");
    GPInstanceRunGP(inst);
    GPDAG& dag = inst.GetDAG();
    GPEngine& gpengine = inst.GetGPEngine();
    inst.MakeNNIEngine();
    NNIEngine& nni_engine = inst.GetNNIEngine();
    // Run unmodified DAG with resized GPEngine test.
    if (perform_resize_unmodded_test) {
      // Verify engine not resized yet by accessing too big index.
      PVId plv_idx_out_of_range =
          PVId((dag.NodeCountWithoutDAGRoot() * 10 * PLVNodeHandler::plv_count_) - 1);
      CHECK_THROWS(gpengine.GetPLV(plv_idx_out_of_range));
      // Force bigger reallocation, with no reindexing.
      gpengine.GrowPLVs(pre_dag.NodeCountWithoutDAGRoot(), std::nullopt,
                        pre_dag.NodeCountWithoutDAGRoot() * 10);
      gpengine.GrowGPCSPs(pre_dag.EdgeCountWithLeafSubsplits(), std::nullopt,
                          pre_dag.EdgeCountWithLeafSubsplits() * 10);
      // Verify engine was resized by accessing too big index.
      CHECK_NOTHROW(gpengine.GetPLV(plv_idx_out_of_range));
      bool gp_run_passes = CheckGPEngineRun(inst, pre_inst);
      return gp_run_passes;
    }
    // Initialize Maps and Reindexers.
    Reindexer node_reindexer, node_reindexer_without_root, edge_reindexer;
    node_reindexer = Reindexer::IdentityReindexer(inst.GetDAG().NodeCount());
    edge_reindexer =
        Reindexer::IdentityReindexer(inst.GetDAG().EdgeCountWithLeafSubsplits());
    // Add NNIs to DAG and check resized and reindexed properly.
    nni_engine.SyncAdjacentNNIsWithDAG();
    size_t nni_count = nni_engine.GetAdjacentNNICount();
    size_t nni_add = 0;
    while (nni_count > 0) {
      for (size_t i = 0; i < nni_count; i++) {
        auto nni = GetWhichNNIFromSet(nni_engine.GetAdjacentNNIs(), i);
        auto mods = inst.GetDAG().AddNodePair(nni.parent_, nni.child_);
        node_reindexer = node_reindexer.ComposeWith(mods.node_reindexer);
        edge_reindexer = edge_reindexer.ComposeWith(mods.edge_reindexer);
        nni_add++;
        if (nni_add >= nni_add_limit) {
          break;
        }
        if (nni_add % test_after_every == 0) {
          node_reindexer_without_root =
              node_reindexer.RemoveNewIndex(dag.GetDAGRootNodeId().value_);
          size_t node_count = dag.NodeCountWithoutDAGRoot();
          size_t edge_count = dag.EdgeCountWithLeafSubsplits();
          if (!skip_reindexing) {
            gpengine.GrowPLVs(node_count, node_reindexer_without_root);
            gpengine.GrowGPCSPs(edge_count, edge_reindexer);
          } else {
            gpengine.GrowPLVs(node_count);
            gpengine.GrowGPCSPs(edge_count);
          }

          // Test resizing and reindexing.
          test_passes =
              CheckGPEngineResizeAndReindex(dag, gpengine, pre_dag, pre_gpengine);
          test_array.push_back(test_passes);
          // Reinitialize reindexers.
          node_reindexer = Reindexer::IdentityReindexer(dag.NodeCount());
          edge_reindexer =
              Reindexer::IdentityReindexer(dag.EdgeCountWithLeafSubsplits());
        }
      }
      nni_engine.ResetAllNNIs();
      nni_engine.SyncAdjacentNNIsWithDAG();
      nni_count = nni_engine.GetAdjacentNNICount();

      if (nni_add >= nni_add_limit) {
        break;
      }
    }
    // Test final resizing and reindexing.
    node_reindexer_without_root =
        node_reindexer.RemoveNewIndex(dag.GetDAGRootNodeId().value_);
    if (!skip_reindexing) {
      gpengine.GrowPLVs(dag.NodeCountWithoutDAGRoot(), node_reindexer_without_root);
      gpengine.GrowGPCSPs(dag.EdgeCountWithLeafSubsplits(), edge_reindexer);
    } else {
      gpengine.GrowPLVs(dag.NodeCountWithoutDAGRoot());
      gpengine.GrowGPCSPs(dag.EdgeCountWithLeafSubsplits());
    }
    test_passes = CheckGPEngineResizeAndReindex(dag, gpengine, pre_dag, pre_gpengine);
    test_array.push_back(test_passes);
    // Finally, test run full GP Optimization after all modifications completed.
    inst.ReinitializePriors();
    dag.ReinitializeTidyVectors();
    inst.EstimateBranchLengths(0.0001, 100, true);
    inst.PopulatePLVs();
    inst.ComputeLikelihoods();
    inst.ComputeMarginalLikelihood();

    test_passes = std::accumulate(test_array.begin(), test_array.end(), true,
                                  std::logical_and<>());
    return test_passes;
  };

  // TEST_0: Test that resize and reindex GPEngine works with no modification the DAG.
  auto test_0 = ResizeAndReindexGPEngineTest(0, 1, false, false);
  CHECK_MESSAGE(test_0,
                "TEST_0: Resize and reindex GPEngine fails when no modifications are "
                "made to DAG.");
  // TEST_1: Test resize and reindex GPEngine works when adding a single node pair to
  // DAG.
  auto test_1 = ResizeAndReindexGPEngineTest(1, 1, false, false);
  CHECK_MESSAGE(
      test_1,
      "TEST_1: Resize and reindex GPEngine fails after single AddNodePair to DAG.");
  // TEST_2: Test that improper mapping occurs when not reindexing GPEngine when adding
  // a single node pair to DAG.
  auto test_2 = ResizeAndReindexGPEngineTest(10, 1, true, false);
  CHECK_FALSE_MESSAGE(test_2,
                      "TEST_2: Resize and reindex GPEngine is not incorrect when NOT "
                      "reindexing after single AddNodePair to DAG.");
  // TEST_3: Test resize and reindex GPEngine works when adding a many node pairs,
  // performing resizing and reindexing for each modification of DAG.
  auto test_3 = ResizeAndReindexGPEngineTest(100, 1, false, false);
  CHECK_MESSAGE(test_3,
                "TEST_3: Resize and reindex GPEngine fails after multiple AddNodePair, "
                "reindexed individually.");
  // TEST_4: Test resize and reindex GPEngine works when adding a many node pairs,
  // composing multiple modifications of DAG into single reindexing operation.
  auto test_4 = ResizeAndReindexGPEngineTest(100, 10, false, false);
  CHECK_MESSAGE(
      test_4,
      "TEST_4: Resize and reindex GPEngine fails after multiple AddNodePair to DAG, "
      "reindexed in batches.");

  // TEST_5: Resizes GPEngine without modifying the DAG.  Then tests that resized
  // GPEngine and unmodified GPEngine produce same GP run results.
  auto test_5 = ResizeAndReindexGPEngineTest(1, 1, true, true);
  CHECK_MESSAGE(
      test_5,
      "TEST_5: Resized GPEngine with unmodified DAG changed results of GP Run.");
};

// Compares NNI likelihoods computed by two different GPInstances.
// - The true GPInstance adds each NNI individually to the DAG, resizes the GPEngine,
// repopulates all PLVs in DAG, then recomputes the likelihoods for each NNI.
// - The NNIEngine version of the GPInstance adds all NNIs to the GraftDAG, then resizes
// GPEngine for temporary space (two PLVs total, plus one BranchLength and
// PerGPCSPLogLikelihood per NNI), then generates and processes a GPOperationVector for
// computing all NNI Likelihoods in series.
// - Results of each version are then compared by the likelihood of each NNI's central
// edge, spanning the added parent and child.  All priors are set to 1.0 and all branch
// lengths set to 0.1 to remove their impact on the computation.
// - Note: Input DAG is fully connected -- all legal edges between any two subsplits in
// DAG are added.  This ensures that NNIs via truthDAG using AddNodePair and graftDAG
// using Pre-NNI references have the same topology.
TEST_CASE("NNIEngine via GPEngine: Proposed NNI vs DAG NNI GPLikelihoods") {
  // Fetch likelihood from instance.
  auto GPInstGetNNILikelihood = [](GPInstance& inst, const NNIOperation& nni) {
    const GPDAG& dag = inst.GetDAG();
    const auto edge_idx = dag.GetEdgeIdx(nni.parent_, nni.child_);
    Assert(edge_idx < size_t(inst.GetGPEngine().GetPerGPCSPLogLikelihoods().size()),
           "Edge idx out of range for GPInstGetNNILikelihood.");
    const double likelihood =
        inst.GetGPEngine().GetPerGPCSPLogLikelihoods()[edge_idx.value_];
    return likelihood;
  };

  auto CompareGPLikelihoodsOfProposedNNIsVsDAGNNIs =
      [&GPInstGetNNILikelihood](const std::string& fasta_path,
                                const std::string& newick_path,
                                const bool do_optimize_new_branch_lengths = false,
                                const bool do_fixed_branch_lengths = false) {
        bool passes_test = true;
        bool do_optimize_branch_lengths = !do_fixed_branch_lengths;

        // Likelihoods.
        NNIScoreMap prenni_predag_likelihoods;
        NNIScoreMap prenni_graftdag_likelihoods;
        NNIScoreMap nni_graftdag_likelihoods;
        NNIScoreMap prenni_truthdag_likelihoods;
        NNIScoreMap nni_truthdag_likelihoods;
        // Branch Lengths.
        BranchLengths predag_branchlengths;
        NNIBranchLengthsMap graftdag_branchlengths;
        NNIBranchLengthsMap truthdag_branchlengths;
        // Branch Map.
        BranchMap predag_branchmap;
        NNIBranchMapMap graftdag_branchmaps;
        NNIBranchMapMap truthdag_branchmaps;

        // PreDAG Instance: unaltered DAG.
        auto pre_inst =
            GPInstanceOfFiles(fasta_path, newick_path, "_ignore/mmapped_pv_pre.data");
        pre_inst.MakeNNIEngine();
        GPDAG& pre_dag = pre_inst.GetDAG();
        GPEngine& pre_gpengine = pre_inst.GetGPEngine();
        NNIEngine& nniengine = pre_inst.GetNNIEngine();
        pre_dag.FullyConnect();
        pre_gpengine.GrowPLVs(pre_dag.NodeCountWithoutDAGRoot());
        pre_gpengine.GrowGPCSPs(pre_dag.EdgeCountWithLeafSubsplits());
        pre_gpengine.SetNullPrior();
        const bool do_init_branch_lengths = false;
        GPInstanceRunGP(pre_inst,
                        (do_optimize_branch_lengths || do_init_branch_lengths),
                        do_fixed_branch_lengths, false);

        predag_branchmap =
            pre_gpengine.GetBranchLengthHandler().BuildBranchLengthMap(pre_dag);
        const auto pre_branches =
            pre_gpengine.GetBranchLengths(0, pre_gpengine.GetPaddedGPCSPCount());
        predag_branchlengths = pre_branches;
        nniengine.SyncAdjacentNNIsWithDAG();

        // Map from pre-NNI to NNI that created NNI.
        std::map<NNIOperation, NNIOperation> nni_to_prenni_map;
        for (const auto& nni : nniengine.GetAdjacentNNIs()) {
          auto pre_nni = nniengine.GetDAG().FindNNINeighborInDAG(nni);
          nni_to_prenni_map.insert({nni, pre_nni});
        }

        // Compute likelihoods for preDAG.
        for (const auto& [nni, pre_nni] : nni_to_prenni_map) {
          std::ignore = nni;
          const auto likelihood = GPInstGetNNILikelihood(pre_inst, pre_nni);
          prenni_predag_likelihoods.insert({pre_nni, likelihood});
        }

        // Instance that is used by graftDAG.
        auto graft_inst =
            GPInstanceOfFiles(fasta_path, newick_path, "_ignore/mmapped_pv_graft.data");
        graft_inst.MakeNNIEngine();
        GPDAG& graft_dag = graft_inst.GetDAG();
        GPEngine& graft_gpengine = graft_inst.GetGPEngine();
        NNIEngine& graft_nniengine = graft_inst.GetNNIEngine();
        graft_dag.FullyConnect();
        graft_gpengine.GrowPLVs(graft_dag.NodeCountWithoutDAGRoot());
        graft_gpengine.GrowGPCSPs(graft_dag.EdgeCountWithLeafSubsplits());
        graft_gpengine.SetNullPrior();
        graft_inst.GetNNIEngine().GetGPEvalEngine().SetOptimizeNewEdges(
            do_optimize_new_branch_lengths);
        graft_gpengine.GetBranchLengthHandler().ApplyBranchLengthMap(predag_branchmap,
                                                                     graft_dag);
        GPInstanceRunGP(graft_inst, false, do_fixed_branch_lengths, false);

        // Compute likelihoods for graftDAG.
        graft_nniengine.SyncAdjacentNNIsWithDAG();
        graft_nniengine.GraftAdjacentNNIsToDAG();
        graft_nniengine.GrowEvalEngineForAdjacentNNIs(true, true);
        graft_gpengine.SetNullPrior();
        graft_nniengine.GetGPEvalEngine().SetOptimizationMaxIteration(20);
        for (const auto& [nni, pre_nni] : nni_to_prenni_map) {
          std::ignore = nni;
          auto pre_nni_llh = GPInstGetNNILikelihood(graft_inst, pre_nni);
          prenni_graftdag_likelihoods.insert({pre_nni, pre_nni_llh});
        }

        for (const auto& [nni, pre_nni] : nni_to_prenni_map) {
          // GraftDAG Instance.
          const auto [graft_llh, _] =
              graft_nniengine.GetGPEvalEngine().ComputeAdjacentNNILikelihood(nni);
          std::ignore = _;
          nni_graftdag_likelihoods.insert({nni, graft_llh});
          const auto graft_branchmap =
              graft_gpengine.GetBranchLengthHandler().BuildBranchLengthMap(graft_dag);
          graftdag_branchmaps[nni] = graft_branchmap;
          const auto graft_branches =
              graft_gpengine.GetBranchLengths(0, graft_gpengine.GetPaddedGPCSPCount());
          graftdag_branchlengths[nni] = graft_branches;

          // TruthDAG Instance.
          auto truth_inst = GPInstanceOfFiles(fasta_path, newick_path,
                                              "_ignore/mmapped_pv_truth.data");
          truth_inst.MakeNNIEngine();
          auto& truth_dag = truth_inst.GetDAG();
          auto& truth_nniengine = truth_inst.GetNNIEngine();
          truth_dag.FullyConnect();
          auto& truth_gpengine = truth_inst.GetGPEngine();
          auto mods = truth_dag.AddNodePair(nni);
          auto node_reindexer_without_root =
              mods.node_reindexer.RemoveNewIndex(truth_dag.GetDAGRootNodeId().value_);
          truth_gpengine.GrowPLVs(truth_dag.NodeCountWithoutDAGRoot(),
                                  node_reindexer_without_root);
          truth_gpengine.GrowGPCSPs(truth_dag.EdgeCountWithLeafSubsplits(),
                                    mods.edge_reindexer);
          truth_gpengine.SetNullPrior();
          GPInstanceRunGP(truth_inst, false, do_fixed_branch_lengths, false);
          truth_gpengine.GetBranchLengthHandler().ApplyBranchLengthMap(predag_branchmap,
                                                                       truth_dag);
          if (!do_fixed_branch_lengths) {
            truth_nniengine.GetGPEvalEngine().CopyGPEngineDataAfterAddingNNI(pre_nni,
                                                                             nni);
          }
          if (do_optimize_new_branch_lengths) {
            truth_nniengine.GetGPEvalEngine().NNIBranchLengthOptimization(nni);
          }
          truth_inst.PopulatePLVs();
          truth_inst.ComputeLikelihoods();
          truth_inst.ComputeMarginalLikelihood();

          const auto truth_branchmap =
              truth_gpengine.GetBranchLengthHandler().BuildBranchLengthMap(truth_dag);
          truthdag_branchmaps[nni] = truth_branchmap;
          const auto truth_branches =
              truth_gpengine.GetBranchLengths(0, truth_gpengine.GetGPCSPCount());
          truthdag_branchlengths[nni] = truth_branches;

          auto pre_nni_llh = GPInstGetNNILikelihood(truth_inst, pre_nni);
          prenni_truthdag_likelihoods.insert({pre_nni, pre_nni_llh});
          auto truth_llh = GPInstGetNNILikelihood(truth_inst, nni);
          nni_truthdag_likelihoods.insert({nni, truth_llh});
        }

        // Tests that pre-NNIs that created new NNIs were unaltered.
        const double tolerance = 1e-3;
        for (const auto& [nni, pre_nni] : nni_to_prenni_map) {
          std::ignore = nni;
          const auto prenni_truth = prenni_truthdag_likelihoods.at(pre_nni);
          const auto prenni_graft = prenni_graftdag_likelihoods.at(pre_nni);
          const auto diff = std::abs(prenni_truth - prenni_graft);
          const auto passes_current_test = (diff < tolerance);
          passes_test = (passes_test & passes_current_test);
          CHECK_MESSAGE(diff < tolerance,
                        "Pre-NNI Likelihood from NNI Engine does not match truth.");
        }
        // Tests that adding new NNIs via GraftDAG produces same likelihood as
        // TruthDAG.
        for (const auto& [nni, pre_nni] : nni_to_prenni_map) {
          std::ignore = pre_nni;
          const auto nni_truth = nni_truthdag_likelihoods.at(nni);
          const auto nni_graft = nni_graftdag_likelihoods.at(nni);
          const auto diff = std::abs(nni_truth - nni_graft);
          const auto passes_current_test = (diff < tolerance);
          if (!passes_current_test) {
            std::cout << "FAIL " << nni_truth << " " << nni_graft << std::endl;
          }
          passes_test = (passes_test & passes_current_test);
          CHECK_MESSAGE(diff < tolerance,
                        "NNI Likelihood from NNI Engine does not match truth.");
        }
        return passes_test;
      };

  // Test_0
  const std::string fasta_path_0 = "data/hello.fasta";
  const std::string newick_path_0 = "data/hello_rooted_diff_branches.nwk";
  CHECK_MESSAGE(CompareGPLikelihoodsOfProposedNNIsVsDAGNNIs(fasta_path_0, newick_path_0,
                                                            false, true),
                "Test_0a: Hello (with default branch lengths) failed.");
  CHECK_MESSAGE(CompareGPLikelihoodsOfProposedNNIsVsDAGNNIs(fasta_path_0, newick_path_0,
                                                            false, false),
                "Test_0b: Hello (without optimized branch lengths) failed.");
  CHECK_MESSAGE(CompareGPLikelihoodsOfProposedNNIsVsDAGNNIs(fasta_path_0, newick_path_0,
                                                            true, false),
                "Test_0c: Hello (with optimized branch lengths) failed.");

  // Test_1
  const std::string fasta_path_1 = "data/six_taxon_longer.fasta";
  const std::string newick_path_1 = "data/six_taxon_rooted_simple.nwk";
  CHECK_MESSAGE(CompareGPLikelihoodsOfProposedNNIsVsDAGNNIs(fasta_path_1, newick_path_1,
                                                            false, true),
                "Test_1a: Six Taxon (with default branch lengths) failed.");
  CHECK_MESSAGE(CompareGPLikelihoodsOfProposedNNIsVsDAGNNIs(fasta_path_1, newick_path_1,
                                                            false, false),
                "Test_1b: Six Taxon (without optimized branch lengths) failed.");
  CHECK_MESSAGE(CompareGPLikelihoodsOfProposedNNIsVsDAGNNIs(fasta_path_1, newick_path_1,
                                                            false, false),
                "Test_1c: Six Taxon (with optimized branch lengths) failed.");
}

// This checks that potential parent or child nodes found via subsplit map lookup match
// those found via brute force (linear scan) after adding or grafting nodes to DAG.
// Repeats tests after adding or grafting NNIs to DAG.
TEST_CASE("NNIEngine: Finding Parent and Child Nodes After Adding/Grafting Nodes") {
  bool is_quiet = false;
  std::stringstream dev_null;
  std::ostream& os = (is_quiet ? dev_null : std::cerr);

  const std::string fasta_path = "data/five_taxon.fasta";
  const std::string newick_path_1 = "data/five_taxon_rooted.nwk";
  auto inst_1 = GPInstanceOfFiles(fasta_path, newick_path_1, "_ignore/mmapped_pv.data");
  inst_1.MakeTPEngine();
  inst_1.MakeNNIEngine();
  auto& dag_1 = inst_1.GetDAG();
  auto& nniengine_1 = inst_1.GetNNIEngine();
  auto& graftdag_1 = nniengine_1.GetGraftDAG();
  nniengine_1.SetTPLikelihoodCutoffFilteringScheme(0.0);
  nniengine_1.SetTopNScoreFilteringScheme(2);
  // nniengine_1.SetNoFilter();
  nniengine_1.RunInit();

  auto VectorToSet = [](const std::pair<NodeIdVector, NodeIdVector>& node_id_vector)
      -> std::pair<std::set<NodeId>, std::set<NodeId>> {
    auto [left, right] = node_id_vector;
    std::set<NodeId> left_set(left.begin(), left.end());
    std::set<NodeId> right_set(right.begin(), right.end());
    return {left_set, right_set};
  };

  auto TestDAGFindNodeIdsViaMapVSViaScan = [&os, &dag_1,
                                            &VectorToSet](const Bitset& subsplit) {
    bool test_passed = true;
    {
      auto [left, right] = VectorToSet(dag_1.FindParentNodeIds(subsplit));
      auto [left_map, right_map] = VectorToSet(dag_1.FindParentNodeIdsViaMap(subsplit));
      auto [left_scan, right_scan] =
          VectorToSet(dag_1.FindParentNodeIdsViaScan(subsplit));
      bool results_match = (left == left_scan and left_map == left_scan) and
                           (right == right_scan and right_map == right_scan);
      if (!results_match) {
        os << "DAG_FIND_PARENT_NODES: " << (results_match ? "PASS" : "FAIL")
           << std::endl;
        os << "LEFT: " << left_map << " " << left_scan << std::endl;
        os << "RIGHT: " << right_map << " " << right_scan << std::endl;
      }
      test_passed &= results_match;
    }
    {
      auto [left, right] = VectorToSet(dag_1.FindChildNodeIds(subsplit));
      auto [left_map, right_map] = VectorToSet(dag_1.FindChildNodeIdsViaMap(subsplit));
      auto [left_scan, right_scan] =
          VectorToSet(dag_1.FindChildNodeIdsViaScan(subsplit));
      bool results_match = (left == left_scan and left_map == left_scan) and
                           (right == right_scan and right_map == right_scan);
      if (!results_match) {
        os << "DAG_FIND_CHILD_NODES: " << (results_match ? "PASS" : "FAIL")
           << std::endl;
        os << "LEFT: " << left_map << " " << left_scan << std::endl;
        os << "RIGHT: " << right_map << " " << right_scan << std::endl;
      }
      test_passed &= results_match;
    }
    return test_passed;
  };

  auto TestGraftDAGFindNodeIdsViaMapVSViaScan = [&os, &graftdag_1,
                                                 &VectorToSet](const Bitset& subsplit) {
    bool test_passed = true;
    {
      auto [left, right] = VectorToSet(graftdag_1.FindParentNodeIds(subsplit));
      auto [left_map, right_map] =
          VectorToSet(graftdag_1.FindParentNodeIdsViaMap(subsplit));
      auto [left_scan, right_scan] =
          VectorToSet(graftdag_1.FindParentNodeIdsViaScan(subsplit));
      bool results_match = (left == left_scan and left_map == left_scan) and
                           (right == right_scan and right_map == right_scan);
      if (!results_match) {
        os << "GRAFTDAG_FIND_PARENT_NODES: " << (results_match ? "PASS" : "FAIL")
           << std::endl;
        os << "LEFT: " << left_map << " " << left_scan << std::endl;
        os << "RIGHT: " << right_map << " " << right_scan << std::endl;
      }
      test_passed &= results_match;
    }
    {
      auto [left, right] = VectorToSet(graftdag_1.FindChildNodeIds(subsplit));
      auto [left_map, right_map] =
          VectorToSet(graftdag_1.FindChildNodeIdsViaMap(subsplit));
      auto [left_scan, right_scan] =
          VectorToSet(graftdag_1.FindChildNodeIdsViaScan(subsplit));
      bool results_match = (left == left_scan and left_map == left_scan) and
                           (right == right_scan and right_map == right_scan);
      if (!results_match) {
        os << "GRAFTDAG_FIND_CHILD_NODES: " << (results_match ? "PASS" : "FAIL")
           << std::endl;
        os << "LEFT: " << left_map << " " << left_scan << std::endl;
        os << "RIGHT: " << right_map << " " << right_scan << std::endl;
      }
      test_passed &= results_match;
    }
    return test_passed;
  };

  for (NodeId node_id{0}; node_id < dag_1.NodeCount(); node_id++) {
    auto subsplit = dag_1.GetDAGNodeBitset(node_id);
    CHECK_MESSAGE(TestDAGFindNodeIdsViaMapVSViaScan(subsplit),
                  "Finding nodes via map does not match finding nodes via scan (before "
                  "adding NNIs).");
  }

  size_t max_iter = 10;
  for (size_t iter = 0; iter < max_iter; iter++) {
    std::cout << "DAG: " << dag_1.NodeCount() << " "
              << dag_1.EdgeCountWithLeafSubsplits() << std::endl;
    nniengine_1.GraftAdjacentNNIsToDAG();
    for (NodeId node_id{0}; node_id < graftdag_1.NodeCount(); node_id++) {
      auto subsplit = graftdag_1.GetDAGNodeBitset(node_id);
      CHECK_MESSAGE(
          TestGraftDAGFindNodeIdsViaMapVSViaScan(subsplit),
          "Finding nodes via map does not match finding nodes via scan (before "
          "adding NNIs).");
    }
    nniengine_1.FilterPreUpdate();
    nniengine_1.FilterEvaluateAdjacentNNIs();
    nniengine_1.FilterPostUpdate();
    nniengine_1.FilterProcessAdjacentNNIs();
    nniengine_1.RemoveAllGraftedNNIsFromDAG();
    nniengine_1.AddAcceptedNNIsToDAG();
    for (NodeId node_id{0}; node_id < dag_1.NodeCount(); node_id++) {
      auto subsplit = dag_1.GetDAGNodeBitset(node_id);
      CHECK_MESSAGE(
          TestDAGFindNodeIdsViaMapVSViaScan(subsplit),
          "Finding nodes via map does not match finding nodes via scan (before "
          "adding NNIs).");
    }
    nniengine_1.RunPostLoop();
  }
}

// ** TPEngine tests **

// Makes and returns an SBNInstance. Used for truth test vs TPEngine likelihoods.
RootedSBNInstance MakeRootedSBNInstance(const std::string& newick_path,
                                        const std::string& fasta_path,
                                        PhyloModelSpecification& specification,
                                        const bool init_time_trees = true,
                                        const size_t thread_count = 1) {
  RootedSBNInstance inst("demo_instance");
  inst.ReadNewickFile(newick_path, false);
  inst.ReadFastaFile(fasta_path);
  inst.ProcessLoadedTrees();
  // make engine and set phylo model parameters
  inst.PrepareForPhyloLikelihood(specification, thread_count);
  return inst;
};

// Build GPInstance with TPEngine and NNIEngine.
GPInstance MakeGPInstanceWithTPEngine(const std::string& fasta_path,
                                      const std::string& newick_path,
                                      const std::string& mmap_path) {
  // Make GPInstance and TPEngine.
  auto inst = GPInstanceOfFiles(fasta_path, newick_path, mmap_path);
  auto& dag = inst.GetDAG();
  inst.MakeGPEngine();
  auto edge_indexer = dag.BuildEdgeIndexer();
  inst.MakeTPEngine();
  inst.MakeNNIEngine();
  TPEngine& tpengine = inst.GetTPEngine();
  tpengine.GetLikelihoodEvalEngine()
      .GetDAGBranchHandler()
      .GetBranchLengths()
      .FillWithDefault();
  // Set choice map according to subsplit method or pcsp method.
  const bool use_subsplit_method = true;
  inst.TPEngineSetChoiceMapByTakingFirst(use_subsplit_method);
  tpengine.GetLikelihoodEvalEngine().Initialize();
  tpengine.GetLikelihoodEvalEngine().ComputeScores();
  return inst;
}

// Builds maps from tree_id->tree and edge_id->tree_id. TreeEdgeMap is determined by
// extracting the TPEngine's choicemap's topology for the given edge.
std::pair<TreeIdMap, TreeEdgeMap>
BuildTreeIdMapAndTreeEdgeMapFromGPInstanceAndChoiceMap(
    GPInstance& inst, bool apply_branch_lengths = false) {
  TreeIdMap final_tree_id_map;
  TreeIdMap tree_id_map;
  TreeEdgeMap tree_edge_map;
  const auto& tp_engine = inst.GetTPEngine();
  for (EdgeId edge_id(0); edge_id < inst.GetDAG().EdgeCountWithLeafSubsplits();
       edge_id++) {
    auto tree = tp_engine.GetTopTreeWithEdge(edge_id);
    tree_id_map.insert({TreeId(edge_id.value_), tree});
  }
  // Build tree_edge_map from tree_id_map.
  for (const EdgeId edge_id : inst.GetDAG().LeafwardEdgeTraversalTrace(false)) {
    if (inst.GetDAG().IsEdgeRoot(edge_id)) {
      continue;
    }
    // Get tree topology from TPEngine's choice map.
    auto topology = tp_engine.GetTopTopologyWithEdge(edge_id);
    bool is_found = false;
    for (const auto& [tree_id, tree] : tree_id_map) {
      if (topology == tree.Topology()) {
        tree_edge_map[edge_id] = tree_id;
        if (apply_branch_lengths) {
          auto& branch_handler =
              inst.GetTPEngine().GetLikelihoodEvalEngine().GetDAGBranchHandler();
          auto final_tree = DAGBranchHandler::BuildTreeWithBranchLengthsFromTopology(
              inst.GetDAG(), branch_handler, tree.Topology());
          final_tree_id_map.insert({tree_id, final_tree});
        } else {
          final_tree_id_map.insert({tree_id, tree});
        }
        is_found = true;
        break;
      }
    }
    Assert(is_found, "Could not find TPEngine topology in TreeCollection.");
  }
  return std::make_pair(final_tree_id_map, tree_edge_map);
}

// Build map from an edge in DAG to its TPEngine score, where each edge corresponds to
// "top tree" topology according to TPEngine's choicemap.
EdgeScoreMap BuildEdgeTPScoreMapFromInstance(GPInstance& inst,
                                             const TPEvalEngineType score_method) {
  EdgeScoreMap tp_score_map;
  auto& tpengine = inst.GetTPEngine();
  // Populate tree vector and edge map.
  const auto& [tree_id_map, tree_edge_map] =
      BuildTreeIdMapAndTreeEdgeMapFromGPInstanceAndChoiceMap(inst);
  std::ignore = tree_id_map;
  // Compute likelihoods with TPEngine.
  if (score_method == TPEvalEngineType::LikelihoodEvalEngine) {
    auto& llh_engine = tpengine.GetLikelihoodEvalEngine();
    // llh_engine.Initialize();
    llh_engine.ComputeScores();
    for (const auto& [edge_id, tree_id] : tree_edge_map) {
      std::ignore = tree_id;
      auto score = llh_engine.GetTopTreeScoreWithEdge(edge_id);
      tp_score_map[edge_id] = score;
    }
  }
  // Compute parsimonies with TPEngine.
  else if (score_method == TPEvalEngineType::ParsimonyEvalEngine) {
    auto& parsimony_engine = tpengine.GetParsimonyEvalEngine();
    parsimony_engine.Initialize();
    parsimony_engine.ComputeScores();
    for (const auto& [edge_id, tree_id] : tree_edge_map) {
      std::ignore = tree_id;
      auto score = parsimony_engine.GetTopTreeScoreWithEdge(edge_id);
      tp_score_map[edge_id] = score;
    }
  } else {
    Failwith("Given TPEvalEngineType is not valid.");
  }

  return tp_score_map;
};

// Build map from an edge in the DAG to its TPEngine proposed NNI score, where each edge
// is from the "Post-DAG", a DAG which already contains all proposed NNIs.
EdgeScoreMap BuildProposedEdgeTPScoreMapFromInstance(
    GPInstance& inst, GPInstance& post_inst, const TPEvalEngineType score_method) {
  EdgeScoreMap tp_score_map;
  auto& nni_engine = inst.GetNNIEngine();
  nni_engine.SyncAdjacentNNIsWithDAG();
  auto& tpengine = inst.GetTPEngine();
  auto& pre_dag = inst.GetDAG();
  auto& post_dag = post_inst.GetDAG();
  auto pre_node_map = SubsplitDAG::BuildNodeIdMapBetweenDAGs(pre_dag, post_dag);
  auto pre_edge_map = SubsplitDAG::BuildEdgeIdMapBetweenDAGs(pre_dag, post_dag);
  auto post_edge_map = SubsplitDAG::BuildEdgeIdMapBetweenDAGs(post_dag, pre_dag);

  if (score_method == TPEvalEngineType::LikelihoodEvalEngine) {
    auto& llh_engine = tpengine.GetLikelihoodEvalEngine();
    llh_engine.Initialize();
    llh_engine.ComputeScores();
    auto best_edge_map =
        tpengine.BuildBestEdgeMapOverNNIs(nni_engine.GetAdjacentNNIs());
    for (const auto& nni : nni_engine.GetAdjacentNNIs()) {
      const auto pre_nni = tpengine.FindHighestPriorityNeighborNNIInDAG(nni);
      const auto post_edge_id = post_dag.GetEdgeIdx(nni);
      double score =
          llh_engine.GetTopTreeScoreWithProposedNNI(nni, pre_nni, 0, best_edge_map);
      tp_score_map[post_edge_id] = score;
    }
  }
  if (score_method == TPEvalEngineType::ParsimonyEvalEngine) {
    auto& parsimony_engine = tpengine.GetParsimonyEvalEngine();
    parsimony_engine.Initialize();
    parsimony_engine.ComputeScores();
    for (const auto& nni : nni_engine.GetAdjacentNNIs()) {
      const auto& pre_nni = nni_engine.GetDAG().FindNNINeighborInDAG(nni);
      double score = parsimony_engine.GetTopTreeScoreWithProposedNNI(nni, pre_nni);
      auto edge_id = post_dag.GetEdgeIdx(nni);
      tp_score_map[edge_id] = score;
    }
  }

  return tp_score_map;
}

// Build map from an edge in DAG to its score, where each edge corresponds
// to "top tree" topology according to TPEngine's choicemap, and its score is
// computed using a BEAGLE likelihood engine.
EdgeScoreMap BuildEdgeScoreMapFromInstanceUsingBeagleEngine(
    GPInstance& inst, TreeIdMap& tree_id_map, TreeEdgeMap& tree_edge_map) {
  EdgeScoreMap score_map;
  const std::string newick_path = inst.GetNewickSourcePath();
  const std::string fasta_path = inst.GetFastaSourcePath();
  PhyloModelSpecification simple_spec{"JC69", "constant", "strict"};
  auto rooted_sbn_inst = MakeRootedSBNInstance(newick_path, fasta_path, simple_spec);
  auto& beagle_engine = *rooted_sbn_inst.GetEngine()->GetFirstFatBeagle();
  for (auto& [edge_id, tree_id] : tree_edge_map) {
    auto score = beagle_engine.UnrootedLogLikelihood(tree_id_map.at(tree_id));
    score_map[edge_id] = score;
  }

  return score_map;
}

// Build map from an edge in DAG to its score, where each edge corresponds
// to "top tree" topology according to TPEngine's choicemap, and its score is
// computed using a Sankoff Handler parsimony engine.
EdgeScoreMap BuildEdgeScoreMapFromInstanceUsingSankoffHandler(
    GPInstance& inst, TreeIdMap& tree_id_map, TreeEdgeMap& tree_edge_map) {
  EdgeScoreMap score_map;
  auto site_pattern = inst.MakeSitePattern();
  auto sankoff_engine = SankoffHandler(site_pattern, "_ignore/sankoff_handler.data");
  for (auto& [edge_id, tree_id] : tree_edge_map) {
    sankoff_engine.RunSankoff(tree_id_map.at(tree_id).Topology());
    auto score = sankoff_engine.ParsimonyScore();
    score_map[edge_id] = score;
  }

  return score_map;
}

bool TestCompareEdgeScoreMapToCorrectEdgeScoreMap(EdgeScoreMap& score_map_test,
                                                  GPDAG& dag_test,
                                                  EdgeScoreMap& score_map_correct,
                                                  GPDAG& dag_correct,
                                                  const bool is_quiet = true) {
  bool is_equal = true;
  std::stringstream dev_null;
  std::ostream& os = (is_quiet ? dev_null : std::cerr);
  const double tolerance = 1e-3;
  for (const auto [edge_id_test, score_test] : score_map_test) {
    double min_diff = abs(score_test - score_map_correct[edge_id_test]);
    double min_score_correct = score_map_correct[edge_id_test];
    EdgeId min_edge_id_correct = edge_id_test;
    // Check all edge scores in correct map for corresponding match in test map.
    for (const auto [edge_id_correct, score_correct] : score_map_correct) {
      const auto diff = abs(score_test - score_correct);
      if (min_diff > diff) {
        min_diff = diff;
        min_score_correct = score_correct;
        min_edge_id_correct = edge_id_correct;
      }
    }
    if (min_diff > tolerance) {
      is_equal = false;
      os << ":: NNI_SCORE_FAIL :: Error: " << min_diff << std::endl;
      os << "CORRECT: " << min_edge_id_correct << ": " << min_score_correct
         << std::endl;
      os << "TEST: " << edge_id_test << ": " << score_test << std::endl;
    }
    CHECK_MESSAGE(min_diff <= tolerance,
                  "Score of Proposed NNIs in smaller DAG without NNIs did not match "
                  "score from larger DAG.");
  }
  if (!is_equal) {
    os << "NOT_EQUAL: " << std::endl;
    os << "CORRECT: " << score_map_correct << std::endl;
    os << "TEST: " << score_map_test << std::endl;
  }
  return is_equal;
}

// Compares the TPEngine's top tree scores for each given edge in the DAG.
// Tests likelihoods by comparing to the likelihoods from BEAGLE engine.
// Tests parsimonies by comparing to the parsimonies from Sankoff handler.
bool CheckAllTPEngineScores(GPInstance& inst, const bool test_likelihood = true,
                            const bool test_parsimony = true,
                            const bool test_pvs = false, const bool is_quiet = true) {
  bool test_passes = true;
  std::stringstream dev_null;
  std::ostream& os = (is_quiet ? dev_null : std::cerr);
  const double tolerance = 1e5;

  auto& dag = inst.GetDAG();
  auto& gpengine = inst.GetGPEngine();
  auto& tpengine = inst.GetTPEngine();
  auto site_pattern = inst.MakeSitePattern();
  const auto& [tree_id_map, tree_edge_map] =
      BuildTreeIdMapAndTreeEdgeMapFromGPInstanceAndChoiceMap(inst);

  // Check that scores from TPEngine match the correct scores from the individual
  // trees. Note, if the test only contains a single tree, then this amounts to
  // checking if each edge's likelihood matches that one tree.
  auto TestMatchingScores = [&os, is_quiet, &tree_edge_map, &test_passes, &tolerance](
                                const std::string& test_name,
                                TreeScoreMap& correct_tree_score_map,
                                EdgeScoreMap& tp_score_map) {
    for (const auto& [edge_id, tree_id] : tree_edge_map) {
      const auto tp_score = tp_score_map[edge_id];
      const auto correct_score = correct_tree_score_map[tree_id];
      double min_error = abs(tp_score - correct_score);
      CHECK_LT(min_error, tolerance);
      if (min_error > tolerance) {
        test_passes = false;
        os << "::" << test_name << "_FAILURE:: EdgeId: " << edge_id
           << ", TP_Score: " << tp_score_map[edge_id]
           << ", Correct_Score: " << correct_tree_score_map[tree_id]
           << ", Error: " << min_error << std::endl;
      }
    }
    if (!test_passes) {
      os << "TestMatchingScore: " << tp_score_map.size() << std::endl;
      os << "TP_Scores: " << tp_score_map.size() << " " << tp_score_map << std::endl;
      os << "Correct_Score: " << correct_tree_score_map.size() << " "
         << correct_tree_score_map << std::endl;
    }
  };

  // Run tests comparing TPEngine for computing likelihood vs. true tree
  // likelihood computed via BEAGLE Engine.
  if (test_likelihood) {
    TreeScoreMap correct_tree_likelihood_map;
    EdgeScoreMap tp_likelihood_map;
    // BEAGLE Engine for correct tree likelihoods.
    const std::string newick_path = inst.GetNewickSourcePath();
    const std::string fasta_path = inst.GetFastaSourcePath();
    PhyloModelSpecification simple_spec{"JC69", "constant", "strict"};
    auto rooted_sbn_inst = MakeRootedSBNInstance(newick_path, fasta_path, simple_spec);
    auto& beagle_engine = *rooted_sbn_inst.GetEngine()->GetFirstFatBeagle();
    for (const auto& [tree_id, tree] : tree_id_map) {
      auto correct_likelihood = beagle_engine.UnrootedLogLikelihood(tree);
      correct_tree_likelihood_map[tree_id] = correct_likelihood;
    }
    // Check against likelihoods with TPEngine.
    for (const auto& [edge_id, tree_id] : tree_edge_map) {
      std::ignore = tree_id;
      auto likelihood =
          tpengine.GetLikelihoodEvalEngine().GetTopTreeScoreWithEdge(edge_id);
      tp_likelihood_map[edge_id] = likelihood;
    }
    // Check that scores from TPEngine match the correct scores from the individual
    // trees computed by BEAGLE engine.
    TestMatchingScores(std::string("LIKELIHOODS"), correct_tree_likelihood_map,
                       tp_likelihood_map);
    // Compare GP and TP partial vectors. Note, this test is only relevant with single
    // trees, as GP and TP PVs are only equal in the case of single tree DAGs.
    if (test_likelihood && test_pvs) {
      auto& tp_pvs = tpengine.GetLikelihoodPVs();
      auto& gp_pvs = gpengine.GetPLVHandler();
      if (test_pvs) {
        test_passes = (tp_pvs.GetCount() == gp_pvs.GetCount());
        for (const auto& pv_type : PLVTypeEnum::Iterator()) {
          for (EdgeId edge_id = 0; edge_id < dag.EdgeCountWithLeafSubsplits();
               edge_id++) {
            NodeId child_id = dag.GetDAGEdge(edge_id).GetChild();
            bool is_equal =
                (tp_pvs.GetPV(pv_type, edge_id) == gp_pvs.GetPV(pv_type, child_id));
            if (!is_equal) {
              test_passes = false;
              os << "!!! *** NOT EQUAL ***" << std::endl;
              os << "TP_" << tp_pvs.ToString(pv_type, edge_id);
              os << "GP_" << gp_pvs.ToString(pv_type, child_id);
            }
          }
        }
      }
    }
  }

  // Run tests comparing TPEngine infastructure for computing parsimony over DAG via
  // Sankoff vs. true trees via Sankoff Handler.  (Sankoff Handler tests are already
  // been tested on full trees in these doctests using trees with known parsimonies,
  // so we can trust it to generate correct tree parsimonies for testing the
  // TPEngine.)
  if (test_parsimony) {
    TreeScoreMap correct_tree_parsimony_map;
    EdgeScoreMap tp_parsimony_map;
    // Sankoff Handler for correct tree parsimonies.
    SankoffHandler sankoff_engine(site_pattern, "_ignore/mmapped_pv.sankoff.data");
    for (const auto& [tree_id, tree] : tree_id_map) {
      sankoff_engine.RunSankoff(tree.Topology());
      auto correct_parsimony = sankoff_engine.ParsimonyScore(tree.Topology()->Id());
      correct_tree_parsimony_map[tree_id] = correct_parsimony;
    }
    // Compute parsimonies with TPEngine.
    for (const auto& [edge_id, tree_id] : tree_edge_map) {
      std::ignore = tree_id;
      auto parsimony =
          tpengine.GetParsimonyEvalEngine().GetTopTreeScoreWithEdge(edge_id);
      tp_parsimony_map[edge_id] = parsimony;
    }
    // Check that scores from TPEngine match the correct scores from the individual
    // trees computed by Sankoff Handler.
    TestMatchingScores(std::string("PARSIMONY"), correct_tree_parsimony_map,
                       tp_parsimony_map);
    // Compare GP and TP partial vectors. Note, this test is only relevant with single
    // trees, as GP and TP PVs are only equal in the case of single tree DAGs.
    if (test_pvs) {
      auto& tp_pvs = tpengine.GetParsimonyPVs();
      auto& sankoff_pvs = sankoff_engine.GetPSVHandler();
      for (const auto& pv_type : PSVTypeEnum::Iterator()) {
        for (EdgeId edge_id = 0; edge_id < tp_pvs.GetCount(); edge_id++) {
          NodeId child_id = dag.GetDAGEdge(edge_id).GetChild();
          bool is_equal =
              (tp_pvs.GetPV(pv_type, edge_id) == sankoff_pvs.GetPV(pv_type, child_id));
          if (!is_equal) {
            test_passes = false;
            os << "!!! *** NOT EQUAL ***" << std::endl;
            os << "TP_" << tp_pvs.ToString(pv_type, edge_id);
            os << "SANKOFF_" << sankoff_pvs.ToString(pv_type, child_id);
          }
        }
      }
    }
  }

  return test_passes;
}

// Compare TPEngine's top tree score (likelihood or parsimony) in DAG for each edge to
// the collection of verified true tree scores that come from the collection of trees
// that created the DAG. In the single tree cases, can compare the partial vectors of
// verified method vs TPEngine method. Note: The input newick file does not need to
// contain every possible tree expressible in the DAG.  The input tree collection only
// needs to be ordered in terms of score.  The DAG edges will then each be assigned
// according to the best/first tree containing the given edge.
bool TestTPEngineScoresAndPVs(const std::string& fasta_path,
                              const std::string& newick_path,
                              const bool test_likelihood = true,
                              const bool test_parsimony = true,
                              const bool test_pvs = false, const bool is_quiet = true) {
  // Make GPInstance and TPEngine.
  auto inst = GPInstanceOfFiles(fasta_path, newick_path, "_ignore/mmapped_pv.data");
  inst.MakeDAG();
  inst.MakeGPEngine();
  inst.MakeTPEngine();
  auto& dag = inst.GetDAG();
  auto& gpengine = inst.GetGPEngine();
  auto& tpengine = inst.GetTPEngine();
  auto tree_collection = inst.GenerateCompleteRootedTreeCollection();
  auto edge_indexer = dag.BuildEdgeIndexer();
  // Compute likelihoods with GPEngine.
  inst.EstimateBranchLengths(0.00001, 100, true);
  inst.PopulatePLVs();
  inst.ComputeLikelihoods();
  // Perform computations for TPEngine.
  tpengine.SetBranchLengths(gpengine.GetBranchLengths());
  tpengine.SetChoiceMapByTakingFirst(tree_collection, edge_indexer);
  if (test_likelihood) {
    tpengine.GetLikelihoodEvalEngine().Initialize();
    tpengine.GetLikelihoodEvalEngine().ComputeScores();
  }
  if (test_parsimony) {
    tpengine.GetParsimonyEvalEngine().Initialize();
    tpengine.GetParsimonyEvalEngine().ComputeScores();
  }
  // Test resulting scores against those by computed by individual trees.
  return CheckAllTPEngineScores(inst, test_likelihood, test_parsimony, test_pvs,
                                is_quiet);
};

// Initializes a ChoiceMap for a DAG. Then uses a naive method that picks the first
// listed neighbor for each parent, sister, left and right child. Tests that results
// is a valid selection (all edges have mapped valid edge choices, except for root
// and leaves).
// - Tests that invalid TreeMask are found invalid.
// - Tests that invalid TreeMasks causes Topology to throw exception.
// - Creates TreeMask and Node Topology for each edge in DAG, a list of edge ids
// which represent a embedded tree in the DAG.
//    - Tests that TreeMask is valid state for the DAG.
//    - Tests that resulting Topology exists in the DAG.
TEST_CASE("TPEngine: ChoiceMap") {
  const std::string fasta_path = "data/six_taxon_longer.fasta";
  const std::string newick_path = "data/six_taxon_rooted_simple.nwk";
  auto inst = GPInstanceOfFiles(fasta_path, newick_path, "_ignore/mmapped_pv.data");
  auto dag = inst.GetDAG();

  auto choice_map = TPChoiceMap(dag);
  CHECK_FALSE_MESSAGE(choice_map.SelectionIsValid(),
                      "ChoiceMap selection was incorrectly found valid.");
  choice_map.SelectFirstEdge();
  CHECK_MESSAGE(choice_map.SelectionIsValid(),
                "ChoiceMap selection was found invalid.");

  // Test for fail states for invalid TreeMasks.
  TPChoiceMap::TreeMask tree_mask;
  Node::NodePtr topology;
  NodeIdVector tree_nodes;
  bool quiet_errors = true;
  for (const auto edge_id : tree_mask) {
    tree_nodes.push_back(dag.GetDAGEdge(edge_id).GetParent());
    tree_nodes.push_back(dag.GetDAGEdge(edge_id).GetChild());
  }
  // Empty tree.
  CHECK_FALSE_MESSAGE(choice_map.TreeMaskIsValid(tree_mask, quiet_errors),
                      "TreeMask is incorrectly valid when empty.");
  CHECK_THROWS_MESSAGE(choice_map.ExtractTopology(tree_mask),
                       "Tree is incorrectly valid when empty.");
  // Complete DAG.
  for (EdgeId edge_id = EdgeId(0); edge_id < dag.EdgeCountWithLeafSubsplits();
       edge_id++) {
    tree_mask.insert(edge_id);
  }
  CHECK_FALSE_MESSAGE(choice_map.TreeMaskIsValid(tree_mask, quiet_errors),
                      "TreeMask is incorrectly valid with full DAG.");
  CHECK_THROWS_MESSAGE(choice_map.ExtractTopology(tree_mask),
                       "Tree is incorrectly valid when missing root edge.");
  // Tree missing root node.
  tree_mask = choice_map.ExtractTreeMask(EdgeId(0));
  for (const auto edge_id : tree_mask) {
    if (dag.IsEdgeRoot(edge_id)) {
      tree_mask.erase(tree_mask.find(edge_id));
      break;
    }
  }
  CHECK_FALSE_MESSAGE(choice_map.TreeMaskIsValid(tree_mask, quiet_errors),
                      "TreeMask is incorrectly valid when missing root edge.");
  CHECK_THROWS_MESSAGE(choice_map.ExtractTopology(tree_mask),
                       "Tree is incorrectly valid when missing root edge.");
  // Tree missing leaf node.
  tree_mask = choice_map.ExtractTreeMask(EdgeId(0));
  for (const auto edge_id : tree_mask) {
    if (dag.IsEdgeLeaf(edge_id)) {
      tree_mask.erase(tree_mask.find(edge_id));
      break;
    }
  }
  CHECK_FALSE_MESSAGE(choice_map.TreeMaskIsValid(tree_mask, quiet_errors),
                      "TreeMask is incorrectly valid when missing leaf edge.");
  CHECK_THROWS_MESSAGE(choice_map.ExtractTopology(tree_mask),
                       "Tree is incorrectly valid when missing leaf edge.");
  // Tree missing internal edge.
  tree_mask = choice_map.ExtractTreeMask(EdgeId(0));
  for (const auto edge_id : tree_mask) {
    if (!dag.IsEdgeLeaf(edge_id) && !dag.IsEdgeLeaf(edge_id)) {
      tree_mask.erase(tree_mask.find(edge_id));
      break;
    }
  }
  CHECK_FALSE_MESSAGE(choice_map.TreeMaskIsValid(tree_mask, quiet_errors),
                      "TreeMask is incorrectly valid when missing internal edge.");
  CHECK_THROWS_MESSAGE(choice_map.ExtractTopology(tree_mask),
                       "Tree is incorrectly valid when missing internal edge.");
  // Tree contains extra edge.
  tree_mask = choice_map.ExtractTreeMask(EdgeId(0));
  for (EdgeId edge_id = EdgeId(0); edge_id < dag.EdgeCountWithLeafSubsplits();
       edge_id++) {
    const auto contains_edge = (tree_mask.find(edge_id) != tree_mask.end());
    if (!contains_edge) {
      const auto& parent_id = dag.GetDAGEdge(edge_id).GetParent();
      const auto contains_parent = (std::find(tree_nodes.begin(), tree_nodes.end(),
                                              parent_id) != tree_nodes.end());
      const auto& child_id = dag.GetDAGEdge(edge_id).GetChild();
      const auto contains_child = (std::find(tree_nodes.begin(), tree_nodes.end(),
                                             child_id) != tree_nodes.end());
      if (!contains_parent && !contains_child) {
        tree_mask.insert(edge_id);
        break;
      }
    }
  }
  CHECK_FALSE_MESSAGE(choice_map.TreeMaskIsValid(tree_mask, quiet_errors),
                      "TreeMask is incorrectly valid when containing an extra edge.");
  CHECK_THROWS_MESSAGE(choice_map.ExtractTopology(tree_mask),
                       "Tree is incorrectly valid when containing an extra edge.");
  // Valid topology that exists in the DAG.
  // ((x0,x1),(x2,((x3,x4),x5)));
  topology = Node::Join(
      Node::Join(Node::Leaf(0, 6), Node::Leaf(1, 6)),
      Node::Join(
          Node::Join(Node::Leaf(2, 6), Node::Join(Node::Leaf(3, 6), Node::Leaf(4, 6))),
          Node::Leaf(5, 6)));
  CHECK_MESSAGE(dag.ContainsTopology(topology, quiet_errors),
                "Incorrectly could not find topology that exists in DAG.");
  // Valid topology that does not exist in the DAG.
  // (((x0,x1),(x2,x3)),(x4, x5))
  topology = Node::Join(Node::Join(Node::Join(Node::Leaf(0, 6), Node::Leaf(1, 6)),
                                   Node::Join(Node::Leaf(2, 6), Node::Leaf(3, 6))),
                        Node::Join(Node::Leaf(4, 6), Node::Leaf(5, 6)));
  CHECK_FALSE_MESSAGE(dag.ContainsTopology(topology, quiet_errors),
                      "Incorrectly found topology that does not exist in DAG.");
  // Incomplete topology -- does not terminate at a leaf.
  // ((x0,x1),(x2,x3_4),x5)));
  topology = Node::Join(
      Node::Join(Node::Leaf(0, 6), Node::Leaf(1, 6)),
      Node::Join(Node::Join(Node::Leaf(2, 6), Node::Leaf(34, Bitset("000110"))),
                 Node::Leaf(5, 6)));
  CHECK_FALSE_MESSAGE(dag.ContainsTopology(topology, quiet_errors),
                      "Incorrectly found incomplete topology in DAG.");
  // Incomplete topology -- subtree with missing taxa.
  // (x2,((x3,x4),x5))
  topology = Node::Join(
      Node::Join(Node::Leaf(2, 6), Node::Join(Node::Leaf(3, 6), Node::Leaf(4, 6))),
      Node::Leaf(5, 6));
  CHECK_FALSE_MESSAGE(dag.ContainsTopology(topology, quiet_errors),
                      "Incorrectly found subtree topology in DAG.");

  // Test TreeMasks created from all DAG edges result in valid tree.
  for (EdgeId edge_id = EdgeId(0); edge_id < dag.EdgeCountWithLeafSubsplits();
       edge_id++) {
    const auto tree_mask = choice_map.ExtractTreeMask(edge_id);
    const auto topology = choice_map.ExtractTopology(edge_id);
    CHECK_MESSAGE(tree_mask.find(edge_id) != tree_mask.end(),
                  "TreeMask did not contain given central edge.");
    CHECK_MESSAGE(choice_map.TreeMaskIsValid(tree_mask, quiet_errors),
                  "Edge resulted in an invalid TreeMask.");
    CHECK_MESSAGE(dag.ContainsTopology(topology, quiet_errors),
                  "Edge resulted in an invalid Topology not contained in DAG.");
  }
}

// Initializes the TPEngine choice map and retrieves the top tree for each edge in
// DAG. Then finds all trees contained in the DAG and verifies that each top tree
// produced is a tree from the DAG.
TEST_CASE("TPEngine: Initialize TPEngine and ChoiceMap") {
  const std::string fasta_path = "data/six_taxon.fasta";
  const std::string newick_path = "data/six_taxon_rooted_simple.nwk";

  auto inst = GPInstanceOfFiles(fasta_path, newick_path);
  GPDAG& dag = inst.GetDAG();
  inst.EstimateBranchLengths(0.00001, 100, true);
  auto all_trees = inst.GenerateCompleteRootedTreeCollection();
  SitePattern site_pattern = inst.MakeSitePattern();
  TPEngine tpengine = TPEngine(dag, site_pattern, "_ignore/mmapped_pv.tpl.data",
                               "_ignore/mmapped_pv.tpp.data");
  tpengine.InitializeChoiceMap();

  auto TopologyExistsInTreeCollection = [](const Node::NodePtr tree_topology,
                                           RootedTreeCollection& tree_collection) {
    for (const auto& tree : tree_collection.Trees()) {
      if (tree.Topology() == tree_topology) {
        return true;
      }
    }
    return false;
  };

  for (EdgeId edge_id = 0; edge_id < dag.EdgeCountWithLeafSubsplits(); edge_id++) {
    auto top_tree_topology = tpengine.GetTopTopologyWithEdge(edge_id);
    auto exists = TopologyExistsInTreeCollection(top_tree_topology, all_trees);
    CHECK_MESSAGE(exists, "Top Tree does not exist in DAG.");
  }
}

// Builds a TPEngine instance from a set of input trees. Then populates TPEngine's PVs
// and computes the top tree likelihood for each edge in the DAG.  Compares these
// likelihoods against the tree's likelihood computed using BEAGLE engine.
TEST_CASE("TPEngine: TPEngine Likelihood scores vs BEAGLE Likelihood scores") {
  auto TestTPLikelihoods = [](const std::string& fasta_path,
                              const std::string& newick_path,
                              const bool test_pvs) -> bool {
    return TestTPEngineScoresAndPVs(fasta_path, newick_path, true, false, test_pvs,
                                    false);
  };
  // Input files.
  const std::string fasta_hello = "data/hello_short.fasta";
  const std::string newick_hello = "data/hello_rooted.nwk";
  const std::string fasta_six = "data/six_taxon.fasta";
  const std::string newick_six_single = "data/six_taxon_rooted_single.nwk";
  const std::string newick_six_simple = "data/six_taxon_rooted_simple.nwk";
  // Test cases.
  const auto test_1 = TestTPLikelihoods(fasta_hello, newick_hello, true);
  CHECK_MESSAGE(test_1, "Hello Example Single Tree failed.");
  const auto test_2 = TestTPLikelihoods(fasta_six, newick_six_single, true);
  CHECK_MESSAGE(test_2, "Six Taxa Single Tree failed.");
  const auto test_3 = TestTPLikelihoods(fasta_six, newick_six_simple, false);
  CHECK_MESSAGE(test_3, "Six Taxa Multi Tree failed.");
}

// Builds a TPEngine instance from a set of input trees. Then populates TPEngine's PVs
// and computes the top tree parsimony for each edge in the DAG.  Compares these
// likelihoods against the tree's likelihood computed using the `sankoff handler`.
TEST_CASE("TPEngine: TPEngine Parsimony scores vs SankoffHandler Parsimony scores") {
  auto TestTPParsimonies = [](const std::string& fasta_path,
                              const std::string& newick_path,
                              const bool test_pvs) -> bool {
    return TestTPEngineScoresAndPVs(fasta_path, newick_path, false, true, test_pvs,
                                    false);
  };
  // Input files.
  const std::string fasta_ex = "data/parsimony_leaf_seqs.fasta";
  const std::string newick_ex = "data/parsimony_tree_0_score_75.0.nwk";
  const std::string fasta_hello = "data/hello_short.fasta";
  const std::string newick_hello = "data/hello_rooted.nwk";
  const std::string fasta_six = "data/six_taxon.fasta";
  const std::string newick_six_single = "data/six_taxon_rooted_single.nwk";
  const std::string newick_six_simple = "data/six_taxon_rooted_simple.nwk";
  const std::string fasta_five = "data/five_taxon.fasta";
  const std::string newick_five_more = "data/five_taxon_rooted_more.nwk";
  // Test cases.
  const auto test_0 = TestTPParsimonies(fasta_ex, newick_ex, false);
  CHECK_MESSAGE(test_0, "Parsimony Test Case Tree failed.");
  const auto test_1 = TestTPParsimonies(fasta_hello, newick_hello, false);
  CHECK_MESSAGE(test_1, "Hello Example Single Tree failed.");
  const auto test_2 = TestTPParsimonies(fasta_six, newick_six_single, false);
  CHECK_MESSAGE(test_2, "Six Taxa Tree failed.");
  const auto test_3 = TestTPParsimonies(fasta_six, newick_six_simple, false);
  CHECK_MESSAGE(test_3, "Six Taxa Multi Tree failed.");
  const auto test_4 = TestTPParsimonies(fasta_five, newick_five_more, false);
  CHECK_MESSAGE(test_4, "Five Taxa Many Trees failed.");
}

// Creates an instance of TPEngine for two DAGs: DAG_1, a simple DAG, and DAG_2, a DAG
// formed from DAG_1 plus all of its adjacent NNIs. Both DAGs PVs are populated and
// their edge TP likelihoods are computed.  Then DAG_1's adjacent proposed NNI
// likelihoods are computed using only PVs from the pre-NNI already contained in
// DAG_1.
// Finally, we compare the results of the proposed NNIs from DAG_1 with the known
// likelihoods of the actual NNIs already contained in DAG_2. This verifies we
// generate the same result from adding NNIs to the DAG and updating as we do from
// using the pre-NNI computation.
TEST_CASE("TPEngine: Proposed NNI vs DAG NNI vs BEAGLE Likelihood") {
  bool print_failures = false;
  // Build NNIEngine from DAG that does not include NNIs. Compute likelihoods.
  auto CompareProposedNNIvsDAGNNIvsBEAGLE =
      [print_failures](const std::string& fasta_path, const std::string& newick_path,
                       const bool optimize_branch_lengths = true,
                       const bool take_first_branch_lengths = true) {
        bool test_passes = true;
        const double tol = 1e-5;
        // Instance for computing proposed scores.
        auto inst_1 = MakeGPInstanceWithTPEngine(fasta_path, newick_path,
                                                 "_ignore/mmapped_pv.1.data");
        // Instance for computing internal DAG scores for comparison.
        auto inst_2 = MakeGPInstanceWithTPEngine(fasta_path, newick_path,
                                                 "_ignore/mmapped_pv.2.data");

        auto& dag_2 = inst_2.GetDAG();
        auto& nni_engine_2 = inst_2.GetNNIEngine();
        auto& tpengine_2 = inst_2.GetTPEngine();

        // Add all NNIs to post-DAG.
        if (take_first_branch_lengths) {
          inst_1.TPEngineSetBranchLengthsByTakingFirst();
          inst_2.TPEngineSetBranchLengthsByTakingFirst();
          inst_1.TPEngineSetChoiceMapByTakingFirst();
          inst_2.TPEngineSetChoiceMapByTakingFirst();
        }
        inst_1.GetTPEngine().GetLikelihoodEvalEngine().SetOptimizeNewEdges(
            optimize_branch_lengths);
        inst_2.GetTPEngine().GetLikelihoodEvalEngine().SetOptimizeNewEdges(
            optimize_branch_lengths);
        inst_1.GetTPEngine().GetLikelihoodEvalEngine().SetOptimizationMaxIteration(1);
        inst_2.GetTPEngine().GetLikelihoodEvalEngine().SetOptimizationMaxIteration(1);

        nni_engine_2.SetNoFilter(true);
        nni_engine_2.RunInit(true);
        CHECK_MESSAGE(tpengine_2.GetChoiceMap().SelectionIsValid(false),
                      "ChoiceMap is not valid before adding NNIs.");
        nni_engine_2.RunMainLoop(true);
        nni_engine_2.RunPostLoop(true);
        CHECK_MESSAGE(tpengine_2.GetChoiceMap().SelectionIsValid(false),
                      "ChoiceMap is not valid after adding NNIs.");

        // Report all NNIs.
        auto& nni_engine = inst_1.GetNNIEngine();
        nni_engine.SyncAdjacentNNIsWithDAG();

        // Likelihoods and Parsimonies of expanded DAG with proposed NNIs.
        auto likelihoods_post = BuildEdgeTPScoreMapFromInstance(
            inst_2, TPEvalEngineType::LikelihoodEvalEngine);
        // Likelihoods and Parsimonies of base DAG.
        auto likelihoods_pre = BuildEdgeTPScoreMapFromInstance(
            inst_1, TPEvalEngineType::LikelihoodEvalEngine);

        // Likelihoods and Parsimonies of DAG's adjacent proposed NNIs.
        auto likelihoods_proposed = BuildProposedEdgeTPScoreMapFromInstance(
            inst_1, inst_2, TPEvalEngineType::LikelihoodEvalEngine);

        // Compute top tree scores using BEAGLE and Sankoff engines.
        auto [tree_id_map, edge_id_map] =
            BuildTreeIdMapAndTreeEdgeMapFromGPInstanceAndChoiceMap(inst_2, true);
        auto beagle = BuildEdgeScoreMapFromInstanceUsingBeagleEngine(
            inst_2, tree_id_map, edge_id_map);

        std::unordered_map<TreeId, std::string> newick_id_map;
        for (const auto& [tree_id, tree] : tree_id_map) {
          newick_id_map[tree_id] = dag_2.TreeToNewickTree(tree);
        }

        // Compare likelihoods by edge_id.
        for (const auto& nni : nni_engine.GetAdjacentNNIs()) {
          auto edge_id = dag_2.GetEdgeIdx(nni);
          auto score_proposed = likelihoods_proposed[edge_id];
          auto score_dag = likelihoods_post[edge_id];
          auto score_tree = beagle[edge_id];
          bool matches_score_dag = (abs(score_proposed - score_dag) < tol);
          bool matches_score_tree = (abs(score_proposed - score_tree) < tol);
          bool matches_dag_tree = (abs(score_dag - score_tree) < tol);

          test_passes &= matches_score_dag;
          test_passes &= matches_score_tree;
          test_passes &= matches_dag_tree;
          if (print_failures) {
            if (!matches_score_dag or !matches_score_tree) {
              std::cout << "SCORE_FAILED_PROPOSED: " << score_tree << " vs "
                        << score_proposed << ": " << abs(score_proposed - score_dag)
                        << std::endl;
              std::cout << "SCORE_FAILED_SCORE: " << score_tree << " vs " << score_dag
                        << ": " << abs(score_dag - score_tree) << std::endl;
            } else {
              std::cout << "SCORE_PASS: " << score_tree << std::endl;
            }
          }
        }

        return test_passes;
      };

  const std::string fasta_path_0 = "data/hello.fasta";
  const std::string newick_path_0 = "data/hello_rooted_diff_branches.nwk";
  CHECK_MESSAGE(
      CompareProposedNNIvsDAGNNIvsBEAGLE(fasta_path_0, newick_path_0, false, false),
      "Test_0a: hello (with fixed branch lengths) failed.");
  CHECK_MESSAGE(
      CompareProposedNNIvsDAGNNIvsBEAGLE(fasta_path_0, newick_path_0, false, true),
      "Test_0b: hello (without optimized branch lengths) failed.");
  CHECK_MESSAGE(
      CompareProposedNNIvsDAGNNIvsBEAGLE(fasta_path_0, newick_path_0, true, true),
      "Test_0c: hello (with optimized branch lengths) failed.");

  const std::string fasta_path_1 = "data/five_taxon.fasta";
  const std::string newick_path_1 = "data/five_taxon_trees_3_4_diff_branches.nwk";
  CHECK_MESSAGE(
      CompareProposedNNIvsDAGNNIvsBEAGLE(fasta_path_1, newick_path_1, false, false),
      "Test_1a: five_taxon_simple (with fixed branch lengths) failed.");
  CHECK_MESSAGE(
      CompareProposedNNIvsDAGNNIvsBEAGLE(fasta_path_1, newick_path_1, false, true),
      "Test_1b: five_taxon_simple (without optimized branch lengths) failed.");
  CHECK_MESSAGE(CompareProposedNNIvsDAGNNIvsBEAGLE(fasta_path_1, newick_path_1, true),
                "Test_1c: five_taxon_simple (with optimized branch lengths) failed.");
}

// Builds TPEngine from single tree DAG, then run branch length optimization.
// Compares results to GPEngine's branch length optimized on the same tree (GP is
// equivalent to traditional likelihood in the single tree case).
TEST_CASE("TPEngine: Branch Length Optimization") {
  GPInstance inst = MakeHelloGPInstance();
  EigenVectorXd seed_branch_lengths(5);
  seed_branch_lengths << 0, 0.22, 0.113, 0.15, 0.1;
  OptimizationMethod method = OptimizationMethod::BrentOptimization;
  bool is_quiet = true;
  bool track_intermediate_values = false;
  // GPEngine
  GPEngine& gp_engine = inst.GetGPEngine();
  inst.GetGPEngine().SetBranchLengths(seed_branch_lengths);
  inst.EstimateBranchLengths(0.0001, 100, is_quiet, track_intermediate_values, method);
  // TPEngine
  inst.MakeTPEngine();
  TPEngine& tp_engine = inst.GetTPEngine();
  tp_engine.GetLikelihoodEvalEngine().GetDAGBranchHandler().SetBranchLengths(
      seed_branch_lengths);
  inst.TPEngineEstimateBranchLengths(0.0001, 100, is_quiet, track_intermediate_values,
                                     method);
  // Capture Results.
  const auto& gpengine_data = gp_engine.GetBranchLengths();
  const auto& tpengine_data = tp_engine.GetBranchLengths();
  // Compare TPEngine results to GPEngine results.
  double tol = 1e-3;
  double max_diff = 0;
  for (size_t i = 0; i < size_t(gpengine_data.size()); i++) {
    if (double diff = abs(tpengine_data[i] - gpengine_data[i]); max_diff < diff) {
      max_diff = diff;
    }
  }

  // Report results.
  if (max_diff > tol) {
    std::cout << "=== TEST FAILED ===" << std::endl;
    std::cout << "max_diff: " << max_diff << std::endl;
    std::cout << "tpengine_lengths:" << tpengine_data << std::endl;
    std::cout << "gpengine_lengths:" << gpengine_data << std::endl;
    std::cout << "tp_likelihood " << tp_engine.GetTopTreeLikelihood(EdgeId(1))
              << std::endl;
    std::cout << "gp_likelihood: " << gp_engine.GetLogMarginalLikelihood() << std::endl;
    auto gp_lengths = gp_engine.GetBranchLengths();
    tp_engine.GetLikelihoodEvalEngine().GetDAGBranchHandler().SetBranchLengths(
        gp_lengths);
    tp_engine.InitializeScores();
    tp_engine.ComputeScores();
    std::cout << "tp_likelihood_using_gp: " << tp_engine.GetTopTreeLikelihoods()
              << std::endl;
  }
  CHECK_LT(max_diff, tol);
}

// Exports Newicks from DAG and TPEngine.  DAG exports a covering set of trees, which
// should contain every node and edge from the DAG in a (unproven) minimum set of trees.
// TPEngine exports an ordered vector of trees which contains a covering set of trees
// that maintains the relative priority of edges according to the TPEngine's choice map.
// Tests this by using the exported trees to build a new DAG and TPEngine, and checks
// that new DAG and TPEngine match the old DAG and TPEngine. Repeats tests after adding
// NNIs to the DAG.
TEST_CASE("TPEngine: Exporting Newicks") {
  const std::string fasta_path = "data/five_taxon.fasta";
  const std::string newick_path_1 = "data/five_taxon_rooted.nwk";
  const std::string newick_path_2 = "data/five_taxon_rooted_shuffled.nwk";

  auto BuildTopTreeMapViaBruteForce = [](GPInstance& inst) {
    std::vector<RootedTree> tree_map;
    for (TreeId tree_id(0); tree_id < inst.GetTPEngine().GetMaxTreeId(); tree_id++) {
      std::vector<RootedTree> temp_trees;
      for (EdgeId edge_id{0}; edge_id < inst.GetDAG().EdgeCountWithLeafSubsplits();
           edge_id++) {
        if (inst.GetTPEngine().GetTreeSource(edge_id) != tree_id) {
          continue;
        }
        auto tree = inst.GetTPEngine().GetTopTreeWithEdge(edge_id);
        bool tree_found = false;
        for (size_t i = 0; i < temp_trees.size(); i++) {
          if (tree == temp_trees[i]) {
            tree_found = true;
            break;
          }
        }
        if (!tree_found) {
          temp_trees.push_back(tree);
        }
      }
      for (const auto& tree : temp_trees) {
        bool tree_found = false;
        for (const auto& old_tree : tree_map) {
          if (tree == old_tree) {
            tree_found = true;
            break;
          }
        }
        if (!tree_found) {
          tree_map.push_back(tree);
        }
      }
    }
    return tree_map;
  };

  // Export a covering newick file from TPEngine, then build new DAG from that file.
  // Compare to the input DAG.
  auto BuildCoveringNewickAndCompareNewDAG = [](GPInstance& inst_1) {
    const std::string temp_newick_path = "_ignore/temp.newick";
    std::ofstream file_out;
    file_out.open(temp_newick_path);
    file_out << inst_1.GetDAG().ToNewickOfCoveringTopologies() << std::endl;
    file_out.close();

    auto inst_2 = GPInstanceOfFiles(inst_1.GetFastaSourcePath(), temp_newick_path,
                                    "_ignore/mmapped_pv.data1");
    bool dags_equal =
        (SubsplitDAG::Compare(inst_1.GetDAG(), inst_2.GetDAG(), false) == 0);
    return dags_equal;
  };

  // Export a top tree newick file from TPEngine, then build new TPEngine from that
  // file. Compare to the input TPEngine.
  auto BuildTopTreeNewickAndCompareNewTPEngine =
      [&BuildTopTreeMapViaBruteForce](GPInstance& inst_1) {
        std::ofstream file_out;
        // Use internal method for building Newick string.
        const std::string temp_newick_path_1 = "_ignore/temp_1.newick";
        std::string newick_1 = inst_1.GetTPEngine().ToNewickOfTopTopologies();
        file_out.open(temp_newick_path_1);
        file_out << newick_1 << std::endl;
        file_out.close();
        // Use brute force method for building Newick string.
        const std::string temp_newick_path_2 = "_ignore/temp_2.newick";
        std::string newick_2;
        file_out.open(temp_newick_path_2);
        const auto tree_map_2 = BuildTopTreeMapViaBruteForce(inst_1);
        for (const auto& tree : tree_map_2) {
          newick_2 += inst_1.GetDAG().TreeToNewickTopology(tree) + '\n';
          file_out << inst_1.GetDAG().TreeToNewickTopology(tree) << std::endl;
        }
        file_out.close();
        bool newicks_equal = (newick_1 == newick_2);
        if (!newicks_equal) {
          std::cerr << "ERROR: Newicks do not match." << std::endl;
          std::cerr << "NEWICK_1: " << std::endl << newick_1 << std::endl;
          std::cerr << "NEWICK_2: " << std::endl << newick_2 << std::endl;
        }

        // Build new TPEngine and check that old and new engines are equal.
        auto inst_2 = GPInstanceOfFiles(inst_1.GetFastaSourcePath(), temp_newick_path_2,
                                        "_ignore/mmapped_pv.data3");
        inst_2.MakeTPEngine();
        inst_2.MakeNNIEngine();
        bool engines_equal =
            (TPEngine::Compare(inst_1.GetTPEngine(), inst_2.GetTPEngine(), false) == 0);
        if (!engines_equal) {
          std::cerr << "ERROR: Engines do not match." << std::endl;
        }
        return (newicks_equal and engines_equal);
      };

  auto inst_1 =
      GPInstanceOfFiles(fasta_path, newick_path_1, "_ignore/mmapped_pv.data1");
  auto inst_2 =
      GPInstanceOfFiles(fasta_path, newick_path_2, "_ignore/mmapped_pv.data2");
  inst_1.MakeTPEngine();
  inst_1.MakeNNIEngine();
  auto& tp_engine = inst_1.GetTPEngine();
  auto& nni_engine = inst_1.GetNNIEngine();
  inst_2.MakeTPEngine();
  inst_2.MakeNNIEngine();

  CHECK_MESSAGE(
      TPEngine::Compare(inst_1.GetTPEngine(), inst_1.GetTPEngine(), false) == 0,
      "TPEngines not equal to self.");
  CHECK_MESSAGE(SubsplitDAG::Compare(inst_1.GetDAG(), inst_2.GetDAG(), true) == 0,
                "DAGs formed from shuffled Newicks not equal.");
  CHECK_MESSAGE(
      TPEngine::Compare(inst_1.GetTPEngine(), inst_2.GetTPEngine(), true) != 0,
      "TPEngines formed from shuffled Newicks are incorrectly equal.");
  CHECK_MESSAGE(BuildCoveringNewickAndCompareNewDAG(inst_1),
                "DAG built from Covering Newick not equal to the DAG that build it "
                "(before adding NNIs).");
  CHECK_MESSAGE(BuildTopTreeNewickAndCompareNewTPEngine(inst_1),
                "Newick and TPEngine built from Top Tree Newick not equal to the "
                "TPEngine that build it (before adding NNIs).");

  tp_engine.GetLikelihoodEvalEngine().SetOptimizationMaxIteration(5);
  tp_engine.GetLikelihoodEvalEngine().BranchLengthOptimization(false);
  nni_engine.SetTPLikelihoodCutoffFilteringScheme(0.0);
  nni_engine.SetTopNScoreFilteringScheme(1);
  nni_engine.SetReevaluateRejectedNNIs(true);
  nni_engine.RunInit(true);

  for (size_t iter = 0; iter < 10; iter++) {
    nni_engine.RunMainLoop(true);

    if (iter == 7) break;
    // Issue #479: this creates an unknown problem on iteration 7.
    // Error occurs during GetLine() in subsplit_dag_storage.hpp:564.
    // Parent-child vertice pair references a line outside range of DAG.
    // (May be an issue with graft addition/removal process?)

    nni_engine.RunPostLoop(true);

    CHECK_MESSAGE(BuildCoveringNewickAndCompareNewDAG(inst_1),
                  "DAG built from Covering Newick not equal to the DAG that build it"
                  " (after adding NNIs).");
    CHECK_MESSAGE(BuildTopTreeNewickAndCompareNewTPEngine(inst_1),
                  "Newicks and TPEngine built from Top Tree Newick not equal to the "
                  "TPEngine that build it (after adding NNIs).");
  }
}

// PVHandler can be reindexed using two methods.  Either via move-copy, where PVs are
// moved so that the PV vector is ordered according to the order of the node_ids they
// represent in the DAG, or via remapping, where a PV map is updated to point at the
// correct PV when given a node_id, which avoids copying.  Test checks that both
// methods have PV with same values.  Repeats tests after adding NNIs to DAG.
TEST_CASE("TPEngine: Resize and Reindex PV Handler") {
  const std::string fasta_path = "data/five_taxon.fasta";
  const std::string newick_path_1 = "data/five_taxon_rooted.nwk";
  auto inst_1 =
      GPInstanceOfFiles(fasta_path, newick_path_1, "_ignore/mmapped_pv.data1");
  auto inst_2 =
      GPInstanceOfFiles(fasta_path, newick_path_1, "_ignore/mmapped_pv.data2");
  // NNI Engine that uses remapping for reindexing PVs.
  inst_1.MakeTPEngine();
  inst_1.MakeNNIEngine();
  auto& tpengine_1 = inst_1.GetTPEngine();
  auto& nniengine_1 = inst_1.GetNNIEngine();
  auto& pvs_1 = tpengine_1.GetLikelihoodEvalEngine().GetPVs();
  pvs_1.SetUseRemapping(false);
  nniengine_1.SetTPLikelihoodCutoffFilteringScheme(0.0);
  nniengine_1.SetTopNScoreFilteringScheme(1);
  nniengine_1.RunInit();
  // NNI Engine that does NOT use remapping for reindexing PVs.
  inst_2.MakeTPEngine();
  inst_2.MakeNNIEngine();
  auto& tpengine_2 = inst_2.GetTPEngine();
  auto& nniengine_2 = inst_2.GetNNIEngine();
  auto& pvs_2 = tpengine_2.GetLikelihoodEvalEngine().GetPVs();
  pvs_2.SetUseRemapping(false);
  nniengine_2.SetTPLikelihoodCutoffFilteringScheme(0.0);
  nniengine_2.SetTopNScoreFilteringScheme(1);
  nniengine_2.RunInit();

  CHECK_MESSAGE(nniengine_1.GetDAG() == nniengine_2.GetDAG(),
                "DAGs do not match (before adding NNIs).");
  auto pvs_match = pvs_1.Compare(pvs_1, pvs_2, false);
  CHECK_MESSAGE(pvs_match, "PVs do not match (before adding NNIs).");

  size_t max_iter = 10;
  for (size_t iter = 0; iter < max_iter; iter++) {
    bool pvs_match;
    nniengine_1.RunMainLoop();
    nniengine_1.RunPostLoop();
    pvs_match = pvs_1.Compare(pvs_1, pvs_2, true);
    CHECK_FALSE_MESSAGE(
        pvs_match, "PVs incorrectly match (after adding NNIs to only one NNIEngine).");

    nniengine_2.RunMainLoop();
    nniengine_2.RunPostLoop();
    pvs_match = pvs_1.Compare(pvs_1, pvs_2, false);
    CHECK_MESSAGE(pvs_match, "PVs do not match (after adding NNIs).");
  }
}

// ** DAGData tests **

// Builds DAGData vectors from the nodes and edges of a DAG. Checks that data is
// resized and reindexed properly after modifying reference DAG.
TEST_CASE("DAGData: Resize and Reindex") {
  std::string fasta_path = "data/five_taxon.fasta";
  std::string newick_path = "data/five_taxon_rooted.nwk";
  std::string mmap_path_1 = "_ignore/mmap.1.data";
  std::string mmap_path_2 = "_ignore/mmap.2.data";

  GPInstance pre_inst =
      MakeGPInstanceWithTPEngine(fasta_path, newick_path, mmap_path_1);
  auto& pre_dag = pre_inst.GetDAG();
  auto& pre_llh_engine =
      pre_inst.GetNNIEngine().GetTPEvalEngine().GetTPEngine().GetLikelihoodEvalEngine();
  pre_llh_engine.Initialize();

  GPInstance inst = MakeGPInstanceWithTPEngine(fasta_path, newick_path, mmap_path_2);
  auto& dag = inst.GetDAG();
  auto& llh_engine =
      inst.GetNNIEngine().GetTPEvalEngine().GetTPEngine().GetLikelihoodEvalEngine();
  llh_engine.Initialize();

  int default_val = -1;
  DAGNodeIntData node_data(dag, default_val);
  DAGEdgeIntData edge_data(dag, default_val);

  // Resize to fit dag.
  size_t spare_count = 10;
  node_data.Resize(dag.NodeCount(), spare_count, std::nullopt, std::nullopt);
  edge_data.Resize(dag.EdgeCountWithLeafSubsplits(), spare_count, std::nullopt,
                   std::nullopt);

  // Check that counts match the size of the DAG.
  CHECK_EQ(node_data.GetCount(), dag.NodeCount());
  CHECK_EQ(edge_data.GetCount(), dag.EdgeCountWithLeafSubsplits());
  // Check that padded size matches spare_count.
  CHECK_EQ(node_data.GetSpareCount(), spare_count);
  CHECK_EQ(edge_data.GetSpareCount(), spare_count);
  // Check that new data is filled with default.
  for (NodeId i = node_data.GetCount(); i < node_data.GetPaddedCount(); i++) {
    CHECK_EQ(node_data(i), default_val);
  }

  // Check that no exceptions thrown while accessing elements in range.
  // Assign each to their index value.
  for (NodeId i = 0; i < node_data.GetPaddedCount(); i++) {
    CHECK_NOTHROW(node_data(NodeId(i)) = i.value_);
  }
  for (EdgeId i = 0; i < edge_data.GetPaddedCount(); i++) {
    CHECK_NOTHROW(edge_data(EdgeId(i)) = i.value_);
  }

  // Grow DAG by adding all adjacent NNIs.
  DAGNodeIntData pre_node_data(node_data);
  DAGEdgeIntData pre_edge_data(edge_data);
  inst.MakeNNIEngine();
  auto& nni_engine = inst.GetNNIEngine();
  nni_engine.SyncAdjacentNNIsWithDAG();
  auto node_reindexer = Reindexer::IdentityReindexer(dag.NodeCount());
  auto edge_reindexer = Reindexer::IdentityReindexer(dag.EdgeCountWithLeafSubsplits());
  for (const auto& nni : nni_engine.GetAdjacentNNIs()) {
    auto mods = dag.AddNodePair(nni);
    node_reindexer = node_reindexer.ComposeWith(mods.node_reindexer);
    edge_reindexer = edge_reindexer.ComposeWith(mods.edge_reindexer);
  }
  // Resize and Reindex the data vectors.
  node_data.Resize(dag, std::nullopt, node_reindexer);
  edge_data.Resize(dag, std::nullopt, edge_reindexer);
  // Grow engine for data.
  llh_engine.GrowEdgeData(edge_reindexer.size(), edge_reindexer);

  // Check that counts match the size of the DAG.
  CHECK_EQ(node_data.GetCount(), dag.NodeCount());
  CHECK_EQ(edge_data.GetCount(), dag.EdgeCountWithLeafSubsplits());
  // Check that padded size matches spare_count.
  CHECK_EQ(node_data.GetSpareCount(), spare_count);
  CHECK_EQ(edge_data.GetSpareCount(), spare_count);

  // Check that node and edge data was reindexed properly.
  const auto& [node_map, edge_map] =
      BuildNodeAndEdgeMapsFromPreDAGToPostDAG(pre_dag, dag);
  for (const auto& [pre_id, post_id] : node_map) {
    CHECK_EQ(node_data(post_id), pre_node_data(pre_id));
  }
  for (const auto& [pre_id, post_id] : edge_map) {
    CHECK_EQ(edge_data(post_id), pre_edge_data(pre_id));
  }
}

// Checks that two identical Newick trees with different orderings yield the same
// tree.
TEST_CASE("GPInstance: Taxon Sorted Tree Collection") {
  for (const bool sort_taxa : {false, true}) {
    GPInstance inst_1("_ignore/mmap.1.data");
    inst_1.ReadFastaFile("data/three_taxon.fasta");
    inst_1.ReadNewickFile("data/three_taxon_1.nwk", sort_taxa);
    GPInstance inst_2("_ignore/mmap.2.data");
    inst_2.ReadFastaFile("data/three_taxon.fasta");
    inst_2.ReadNewickFile("data/three_taxon_2.nwk", sort_taxa);

    auto& trees_1 = inst_1.GetCurrentlyLoadedTrees();
    auto& trees_2 = inst_2.GetCurrentlyLoadedTrees();
    auto& tree_1 = trees_1.GetTree(0);
    auto& tree_2 = trees_2.GetTree(0);
    if (sort_taxa) {
      CHECK_MESSAGE(tree_1 == tree_2, "Trees incorrectly found not equal.");
    } else {
      CHECK_MESSAGE(tree_1 != tree_2, "Trees incorrectly found equal.");
    }
  }
}
