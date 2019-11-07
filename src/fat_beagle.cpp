// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "fat_beagle.hpp"
#include "beagle.hpp"

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

FatBeagle::FatBeagle(const PhyloModel &phylo_model,
                     const SitePattern &site_pattern)
    : phylo_model_(phylo_model), rescaling_(false) {
  beagle_instance_ = beagle::CreateInstance(site_pattern);
  beagle::PrepareBeagleInstance(beagle_instance_, site_pattern);
  beagle::SetJCModel(beagle_instance_);
};

FatBeagle::~FatBeagle() {
  auto finalize_result = beagleFinalizeInstance(beagle_instance_);
  if (finalize_result != 0) {
    std::cout << "beagleFinalizeInstance gave nonzero return value!";
    std::terminate();
  }
}

Tree PrepareTreeForLikelihood(const Tree &tree) {
  if (tree.Children().size() == 3) {
    return tree.Detrifurcate();
  }  // else
  if (tree.Children().size() == 2) {
    return tree;
  }
  // else
  Failwith(
      "Tree likelihood calculations should be done on a tree with a "
      "bifurcation or a trifurcation at the root.");
}

double FatBeagle::LogLikelihood(
    const Tree &in_tree) {
  beagleResetScaleFactors(beagle_instance_, 0);
  auto tree = PrepareTreeForLikelihood(in_tree);
  auto node_count = tree.BranchLengths().size();
  int int_taxon_count = static_cast<int>(tree.LeafCount());
  std::vector<BeagleOperation> operations;
  tree.Topology()->PostOrder(
      [&operations, int_taxon_count, rescaling = rescaling_](const Node *node) {
        if (!node->IsLeaf()) {
          Assert(node->Children().size() == 2,
                 "Tree isn't bifurcating in LogLikelihood.");
          int dest = static_cast<int>(node->Id());
          int child0_id = static_cast<int>(node->Children()[0]->Id());
          int child1_id = static_cast<int>(node->Children()[1]->Id());
          BeagleOperation op = {
              dest,                            // dest
              BEAGLE_OP_NONE, BEAGLE_OP_NONE,  // src and dest scaling
              child0_id,      child0_id,       // src1 and matrix1
              child1_id,      child1_id        // src2 and matrix2
          };
          if (rescaling) {
            // We don't need scaling buffers for the leaves.
            // Index 0 is reserved for accumulating the sum of log scalers.
            // Thus the scaling buffers are indexed by the edge number minus the
            // number of leaves + 1.
            op.destinationScaleWrite = dest - int_taxon_count + 1;
          }
          operations.push_back(op);
        }
      });
  std::vector<int> node_indices(node_count - 1);
  std::iota(node_indices.begin(), node_indices.end(), 0);
  beagleUpdateTransitionMatrices(beagle_instance_,
                                 0,  // eigenIndex
                                 node_indices.data(),
                                 nullptr,  // firstDerivativeIndices
                                 nullptr,  // secondDervativeIndices
                                 tree.BranchLengths().data(),
                                 static_cast<int>(node_count - 1));

  // This is the entry of scaleBuffer in which we store accumulated factors.
  int cumululative_index = rescaling_ ? 0 : BEAGLE_OP_NONE;
  beagleUpdatePartials(beagle_instance_,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       cumululative_index);  // cumulative scale index

  double log_like = 0;
  std::vector<int> root_id = {static_cast<int>(tree.Id())};
  std::vector<int> category_weight_index = {0};
  std::vector<int> state_frequency_index = {0};
  std::vector<int> cumulative_scale_index = {cumululative_index};
  // We're not exacty sure what this argument is for.
  // The beagle docs say: Number of partialsBuffer to integrate (input)
  // In the BEASTs it's hardcoded to 1 and in MrBayes it appears to be for
  // covarion models.
  int mysterious_count = 1;
  beagleCalculateRootLogLikelihoods(
      beagle_instance_, root_id.data(), category_weight_index.data(),
      state_frequency_index.data(), cumulative_scale_index.data(),
      mysterious_count, &log_like);
  return log_like;
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradient(
    const Tree &in_tree) {
  return beagle::BranchGradient(beagle_instance_, in_tree, rescaling_);
}

double FatBeagle::StaticLogLikelihood(
    FatBeagle *fat_beagle,
    const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "Null FatBeagle pointer!");
  return fat_beagle->LogLikelihood(in_tree);
}

std::pair<double, std::vector<double>> FatBeagle::StaticBranchGradient(
    FatBeagle *fat_beagle,
    const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "Null FatBeagle pointer!");
  return fat_beagle->BranchGradient(in_tree);
}
