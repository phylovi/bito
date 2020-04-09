// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_SBN_INSTANCE_HPP_
#define SRC_ROOTED_SBN_INSTANCE_HPP_

#include "libsbn.hpp"

class RootedSBNInstance : public SBNInstance {
 public:
  explicit RootedSBNInstance(const std::string &name) : SBNInstance(name) {}

  // ** Phylogenetic likelihood

  std::vector<double> LogLikelihoods();
  // For each loaded tree, returns a pair of (likelihood, gradient).
  std::vector<std::pair<double, std::vector<double>>> BranchGradients();

  // ** I/O

  void ReadNewickFile(std::string fname);
  void ReadNexusFile(std::string fname);

  RootedTreeCollection rooted_tree_collection_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedSBNInstance: gradients") {
  RootedSBNInstance inst("charlie");
  inst.ReadNexusFile("data/DS1.BEAST.trees");
  inst.ReadFastaFile("data/DS1.fasta");
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  inst.PrepareForPhyloLikelihood(simple_specification, 2);
  // Issue 187: this should get replaced by node height gradients.
  std::cout << inst.BranchGradients() << std::endl;
}

TEST_CASE("RootedSBNInstance: parsing dates") {
  RootedSBNInstance inst("charlie");
  inst.ReadNexusFile("data/test_beast_tree_parsing.nexus");
  std::vector<double> dates;
  for (auto [tag, date] : inst.rooted_tree_collection_.tag_date_map_) {
    dates.push_back(date);
  }
  std::sort(dates.begin(), dates.end());
  CHECK_EQ(dates[0], 0);
  CHECK_EQ(dates.back(), 80.0);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_ROOTED_SBN_INSTANCE_HPP_
