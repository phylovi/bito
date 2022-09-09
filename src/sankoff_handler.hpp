// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Methods to calculate Sankoff on a tree.

#pragma once

#include "eigen_sugar.hpp"
#include "sugar.hpp"
#include "sankoff_matrix.hpp"
#include "site_pattern.hpp"
#include "node.hpp"
#include "driver.hpp"
#include "pv_handler.hpp"
#include "gp_dag.hpp"

// Partial vector for one node across all sites
using SankoffPartial = NucleotidePLV;
// Each SankoffPartial represents calculations for one node
using SankoffPartialVec = std::vector<SankoffPartial>;
// references for SankoffPartials
using SankoffPartialRef = Eigen::Ref<SankoffPartial>;
using SankoffPartialRefVec = std::vector<SankoffPartialRef>;

class SankoffHandler {
 public:
  // DNA assumption
  static constexpr size_t state_count_ = 4;
  static constexpr double big_double_ = static_cast<double>(INT_MAX);

  // Constructors
  SankoffHandler(SitePattern site_pattern, const std::string &mmap_file_path,
                 double resizing_factor = 2.0)
      : mutation_costs_(SankoffMatrix()),
        site_pattern_(std::move(site_pattern)),
        resizing_factor_(resizing_factor),
        psv_handler_(mmap_file_path, 0, site_pattern_.PatternCount(),
                     resizing_factor_) {
    Assert(site_pattern_.SequenceCount() == site_pattern_.TaxonCount(),
           "Error in SankoffHandler constructor 1: Every sequence should be associated "
           "with a node.");
    psv_handler_.SetNodeCount(site_pattern_.TaxonCount());
    psv_handler_.SetAllocatedNodeCount(
        size_t(ceil(double(psv_handler_.GetPaddedNodeCount()) * resizing_factor_)));
    psv_handler_.Resize(site_pattern_.TaxonCount(),
                        psv_handler_.GetAllocatedNodeCount());
  }

  SankoffHandler(CostMatrix cost_matrix, SitePattern site_pattern,
                 const std::string &mmap_file_path, double resizing_factor = 2.0)
      : mutation_costs_(SankoffMatrix(cost_matrix)),
        site_pattern_(std::move(site_pattern)),
        resizing_factor_(resizing_factor),
        psv_handler_(mmap_file_path, 0, site_pattern_.PatternCount(),
                     resizing_factor_) {
    Assert(site_pattern_.SequenceCount() == site_pattern_.TaxonCount(),
           "Error in SankoffHandler constructor 2: Every sequence should be associated "
           "with a node.");
    psv_handler_.SetNodeCount(site_pattern_.TaxonCount());
    psv_handler_.SetAllocatedNodeCount(
        size_t(ceil(double(psv_handler_.GetPaddedNodeCount()) * resizing_factor_)));
    psv_handler_.Resize(site_pattern_.TaxonCount(),
                        psv_handler_.GetAllocatedNodeCount());
  }

  SankoffHandler(SankoffMatrix sankoff_matrix, SitePattern site_pattern,
                 const std::string &mmap_file_path, double resizing_factor = 2.0)
      : mutation_costs_(std::move(sankoff_matrix)),
        site_pattern_(std::move(site_pattern)),
        resizing_factor_(resizing_factor),
        psv_handler_(mmap_file_path, 0, site_pattern_.PatternCount(),
                     resizing_factor_) {
    Assert(site_pattern_.SequenceCount() == site_pattern_.TaxonCount(),
           "Error in SankoffHandler constructor 3: Every sequence should be associated "
           "with a node.");
    psv_handler_.SetNodeCount(site_pattern_.TaxonCount());
    psv_handler_.SetAllocatedNodeCount(
        size_t(ceil(double(psv_handler_.GetPaddedNodeCount()) * resizing_factor_)));
    psv_handler_.Resize(site_pattern_.TaxonCount(),
                        psv_handler_.GetAllocatedNodeCount());
  }

  PSVHandler &GetPSVHandler() { return psv_handler_; }
  SankoffMatrix &GetCostMatrix() { return mutation_costs_; }

  // Resize PVs to fit model.
  void Resize(const size_t new_node_count);

  // Partial Sankoff Vector Handler.
  SankoffPartialVec PartialsAtPattern(PSVType psv_type, size_t pattern_idx) {
    SankoffPartialVec partials_at_pattern(psv_handler_.GetNodeCount());
    for (NodeId node = 0; node < psv_handler_.GetNodeCount(); node++) {
      partials_at_pattern[node.value_] = psv_handler_(psv_type, node).col(pattern_idx);
    }
    return partials_at_pattern;
  }

  // Fill in leaf-values for P partials.
  void GenerateLeafPartials();

  // Sum p-partials for right and left children of node 'node_id'
  // In this case, we get the full p-partial of the given node after all p-partials
  // have been concatenated into one SankoffPartialVector
  EigenVectorXd TotalPPartial(NodeId node_id, size_t site_idx);

  // Calculate the partial for a given parent-child pair
  EigenVectorXd ParentPartial(EigenVectorXd child_partials);

  // Populate rootward parsimony PV for node.
  void PopulateRootwardParsimonyPVForNode(const NodeId parent_id,
                                          const NodeId left_child_id,
                                          const NodeId right_child_id);
  // Populate leafward parsimony PV for node.
  void PopulateLeafwardParsimonyPVForNode(const NodeId parent_id,
                                          const NodeId left_child_id,
                                          const NodeId right_child_id);

  // Calculates left p_partials, right p_partials, and q_partials for all nodes at all
  // sites in tree topology.
  void RunSankoff(Node::NodePtr topology);

  // Calculates parsimony score on given node across all sites.
  double ParsimonyScore(NodeId node_id);

 private:
  SankoffMatrix mutation_costs_;
  SitePattern site_pattern_;
  double resizing_factor_;
  PSVNodeHandler psv_handler_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("SankoffHandler: Tests on single site sequence.") {
  auto fasta_file = "data/hello_single_nucleotide.fasta";
  auto newick_file = "data/hello_rooted.nwk";
  Alignment alignment = Alignment::ReadFasta(fasta_file);
  Driver driver;
  RootedTreeCollection tree_collection =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(newick_file));
  SitePattern site_pattern = SitePattern(alignment, tree_collection.TagTaxonMap());
  Node::NodePtr topology = tree_collection.GetTree(0).Topology();
  size_t taxon_count = site_pattern.TaxonCount();

  // transitions have cost 1 and transversions have cost 2.5
  auto costs = CostMatrix();
  costs << 0., 2.5, 1., 2.5,  //
      2.5, 0., 2.5, 1.,       //
      1., 2.5, 0., 2.5,       //
      2.5, 1., 2.5, 0.;       //

  SankoffHandler sh = SankoffHandler(costs, site_pattern, "_ignore/mmapped_psv.data");

  // testing RunSankoff for one site
  sh.RunSankoff(topology);

  // testing GenerateLeafPartials() method (which is run as first step of RunSankoff)
  SankoffPartialVec leaves_test = sh.PartialsAtPattern(PSVType::PLeft, 0);
  auto big_double = SankoffHandler::big_double_;
  SankoffPartial leaves_correct_pattern_0(4, topology->Id() + 1);
  // column 1 is G(jupiter), column 2 is C(mars), Column 3 is G(saturn)
  leaves_correct_pattern_0 << big_double, big_double, big_double, 0., 0., big_double,
      0., big_double, 0., 0., 0., big_double, 0., 0., 0., big_double, big_double,
      big_double, 0., 0.;
  for (size_t r = 0; r < taxon_count; r++) {
    CHECK(leaves_test[r].isApprox(leaves_correct_pattern_0.col(r)));
  }

  // test parsimony score for RunSankoff
  CHECK_LT(fabs(sh.ParsimonyScore(0) - 2.5), 1e-10);

  // testing 3rd constructor: SankoffHandler(SankoffMatrix, SitePattern) constructor
  SankoffMatrix sm = SankoffMatrix(costs);
  SankoffHandler sh2 = SankoffHandler(sm, site_pattern, "_ignore/mmapped_psv.data");

  // testing ParentPartial()
  auto child_partials = Eigen::Matrix<double, 4, 2>();
  child_partials << 2.5, 3.5, 3.5, 3.5, 2.5, 3.5, 3.5, 4.5;
  auto parent_test = Eigen::Matrix<double, 4, 1>();
  parent_test.setZero();
  for (size_t child = 0; child < 2; child++) {
    parent_test += sh2.ParentPartial(child_partials.col(child));
  }
  auto parent_correct = Eigen::Matrix<double, 4, 1>();
  parent_correct << 6., 7., 6., 8.;
  CHECK(parent_test.isApprox(parent_correct));
}

TEST_CASE("SankoffHandler: Asymmetric cost matrix test on single site sequence.") {
  auto fasta_file = "data/hello_single_nucleotide.fasta";
  auto newick_file = "data/hello_rooted.nwk";
  Alignment alignment = Alignment::ReadFasta(fasta_file);
  Driver driver;
  RootedTreeCollection tree_collection =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(newick_file));
  SitePattern site_pattern = SitePattern(alignment, tree_collection.TagTaxonMap());
  Node::NodePtr topology = tree_collection.GetTree(0).Topology();

  // transitions have cost 1 and transversions have cost 2.5
  auto costs = CostMatrix();
  costs << 0., 2., 3., 4.,  //
      5., 0., 7., 8.,       //
      9., 10., 0., 12.,     //
      13., 14., 15., 0.;    //

  SankoffHandler sh = SankoffHandler(costs, site_pattern, "_ignore/mmapped_psv.data");
  sh.RunSankoff(topology);
  CHECK_LT(fabs(sh.ParsimonyScore(0) - 8.), 1e-10);
}

TEST_CASE("SankoffHandler: Testing sequence gap characters in GenerateLeafPartials()") {
  auto fasta_file = "data/hello.fasta";
  auto newick_file = "data/hello_rooted.nwk";
  Alignment alignment = Alignment::ReadFasta(fasta_file);
  Driver driver;
  RootedTreeCollection tree_collection =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(newick_file));
  SitePattern site_pattern = SitePattern(alignment, tree_collection.TagTaxonMap());
  Node::NodePtr topology = tree_collection.GetTree(0).Topology();
  size_t taxon_count = site_pattern.TaxonCount();

  // test set up for SankoffHandler with default cost matrix
  SankoffHandler default_sh = SankoffHandler(site_pattern, "_ignore/mmapped_psv.data");

  // testing GenerateLeafPartials() method
  default_sh.GenerateLeafPartials();
  auto leaves_test = default_sh.PartialsAtPattern(PSVType::PLeft, 14);
  SankoffPartial leaves_correct_pattern_14(4, topology->Id() + 1);
  auto big_double = SankoffHandler::big_double_;
  // column 1 is G(jupiter), column 2 is -(mars), Column 3 is G(saturn)
  leaves_correct_pattern_14 << big_double, 0., big_double, 0., 0., big_double, 0.,
      big_double, 0., 0., 0., 0., 0., 0., 0., big_double, 0., big_double, 0., 0.;
  for (size_t r = 0; r < taxon_count; r++) {
    CHECK(leaves_test[r].isApprox(leaves_correct_pattern_14.col(r)));
  }
}

TEST_CASE("SankoffHandler: RunSankoff and ParsimonyScore Tests") {
  auto fasta_file = "data/parsimony_leaf_seqs.fasta";
  auto newick_file = "data/parsimony_tree_0_score_75.0.nwk";
  Alignment alignment = Alignment::ReadFasta(fasta_file);
  Driver driver;
  RootedTreeCollection tree_collection =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(newick_file));
  SitePattern site_pattern = SitePattern(alignment, tree_collection.TagTaxonMap());
  Node::NodePtr topology = tree_collection.GetTree(0).Topology();

  // test set up for SankoffHandler with default cost matrix
  SankoffHandler default_sh = SankoffHandler(site_pattern, "_ignore/mmapped_psv.data");

  double parsimony_score_correct = 75.;
  default_sh.RunSankoff(topology);

  for (NodeId node_id = 0; node_id < topology->Id() + 1; node_id++) {
    CHECK_LT(fabs(default_sh.ParsimonyScore(node_id) - parsimony_score_correct), 1e-10);
  }
}

#endif  // DOCTEST_LIBRARY_INCLUDED
