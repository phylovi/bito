#pragma once

#include "../src/sankoff_handler.hpp"

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
