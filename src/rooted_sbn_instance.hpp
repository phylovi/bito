// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_SBN_INSTANCE_HPP_
#define SRC_ROOTED_SBN_INSTANCE_HPP_

#include "sbn_instance.hpp"

class RootedSBNInstance : public SBNInstance {
 public:
  using SBNInstance::SBNInstance;

  size_t TaxonCount() const override { return tree_collection_.TaxonCount(); }
  StringVector TaxonNames() const override { return tree_collection_.TaxonNames(); }
  size_t TreeCount() const override { return tree_collection_.TreeCount(); }
  TagStringMap TagTaxonMap() const override { return tree_collection_.TagTaxonMap(); }
  Node::TopologyCounter TopologyCounter() const override {
    return tree_collection_.TopologyCounter();
  }
  BitsetSizeDict RootsplitCounterOf(
      const Node::TopologyCounter &topologies) const override {
    return RootedSBNMaps::RootsplitCounterOf(topologies);
  }
  PCSSDict PCSSCounterOf(const Node::TopologyCounter &topologies) const override {
    return RootedSBNMaps::PCSSCounterOf(topologies);
  }
  // Use the loaded trees to get the SBN maps, set taxon_names_, and prepare the
  // sbn_parameters_ vector.
  void ProcessLoadedTrees();

  // ** Phylogenetic likelihood

  std::vector<double> LogLikelihoods();
  // For each loaded tree, returns a pair of (likelihood, gradient).
  std::vector<std::pair<double, std::vector<double>>> BranchGradients();

  // ** I/O

  void ReadNewickFile(std::string fname);
  void ReadNexusFile(std::string fname);

  RootedTreeCollection tree_collection_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedSBNInstance: gradients") {
  RootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/fluA.tree");
  inst.ReadFastaFile("data/fluA.fa");
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  inst.PrepareForPhyloLikelihood(simple_specification, 1);
  auto param_map = inst.GetPhyloModelParamBlockMap();
  param_map.at(StrictClockModel::rate_key_).setConstant(0.001);

  auto likelihood = inst.LogLikelihoods();
  double physher_ll = -4777.616349;
  CHECK_LT(fabs(likelihood[0] - physher_ll), 0.0001);

  auto gradients = inst.BranchGradients();
  std::vector<double> physher_gradients = {
      -0.593654, 6.441290,   11.202945, 5.173924,  -0.904631, 2.731402,   3.157131,
      7.082914,  10.305417,  13.988206, 20.709336, 48.897993, 99.164949,  130.205747,
      17.314019, 21.033290,  -1.336335, 12.259822, 22.887291, 27.176564,  47.487426,
      3.637276,  12.955169,  15.315953, 83.254605, -3.806996, 105.385095, 4.874023,
      22.754466, 6.036534,   25.651478, 29.535185, 29.598789, 1.817247,   10.598685,
      76.259248, 56.481423,  10.679778, 6.587179,  3.330556,  -4.622247,  33.417304,
      63.415767, 188.809515, 23.540875, 17.421076, 1.222568,  22.372012,  34.239511,
      3.486115,  4.098873,   13.200954, 19.726890, 96.808738, 4.240029,   7.414585,
      48.871694, 3.488516,   82.969065, 9.009334,  8.032474,  3.981016,   6.543650,
      53.702423, 37.835952,  2.840831,  7.517186,  19.936861};
  for (size_t i = 0; i < gradients[0].second.size(); i++) {
    CHECK_LT(fabs(gradients[0].second[i] - physher_gradients[i]), 0.0001);
  }
  CHECK_LT(fabs(gradients[0].first - physher_ll), 0.0001);
}

TEST_CASE("RootedSBNInstance: parsing dates") {
  RootedSBNInstance inst("charlie");
  inst.ReadNexusFile("data/test_beast_tree_parsing.nexus");
  std::vector<double> dates;
  for (auto [tag, date] : inst.tree_collection_.tag_date_map_) {
    dates.push_back(date);
  }
  std::sort(dates.begin(), dates.end());
  CHECK_EQ(dates[0], 0);
  CHECK_EQ(dates.back(), 80.0);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_ROOTED_SBN_INSTANCE_HPP_
