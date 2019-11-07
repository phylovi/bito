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
  beagleResetScaleFactors(beagle_instance_, 0);
  auto tree = PrepareTreeForLikelihood(in_tree);
  tree.SlideRootPosition();

  size_t node_count = tree.BranchLengths().size();
  int int_node_count = static_cast<int>(node_count);
  int int_taxon_count = static_cast<int>(tree.LeafCount());
  int internal_count = int_taxon_count - 1;
  std::vector<int> node_indices(node_count - 1);
  std::iota(node_indices.begin(), node_indices.end(), 0);
  std::vector<int> gradient_indices(node_count - 1);
  std::iota(gradient_indices.begin(), gradient_indices.end(), node_count);
  std::vector<BeagleOperation> operations;

  int fixed_node_id = static_cast<int>(tree.Topology()->Children()[1]->Id());
  int root_child_id = static_cast<int>(tree.Topology()->Children()[0]->Id());

  // Calculate lower partials
  tree.Topology()->BinaryIdPostOrder(
      [&operations, int_taxon_count, rescaling = rescaling_](
          int node_id, int child0_id, int child1_id) {
        BeagleOperation op = {
            node_id,                         // dest
            BEAGLE_OP_NONE, BEAGLE_OP_NONE,  // src and dest scaling
            child0_id,      child0_id,       // src1 and matrix1
            child1_id,      child1_id        // src2 and matrix2
        };
        if (rescaling) {
          op.destinationScaleWrite = node_id - int_taxon_count + 1;
        }
        operations.push_back(op);
      });

  // Calculate upper partials
  tree.Topology()->TripleIdPreOrderBifurcating(
      [&operations, &root_child_id, &fixed_node_id, rescaling = rescaling_,
       internal_count,
       int_node_count](int node_id, int sister_id, int parent_id) {
        if (node_id != root_child_id && node_id != fixed_node_id) {
          int upper_partial_index;
          int upper_matrix_index;
          if (parent_id == root_child_id) {
            upper_matrix_index = static_cast<int>(root_child_id);
            upper_partial_index = static_cast<int>(fixed_node_id);
          } else if (parent_id == fixed_node_id) {
            upper_matrix_index = static_cast<int>(root_child_id);
            upper_partial_index = static_cast<int>(root_child_id);
          } else {
            upper_partial_index = parent_id + int_node_count;
            upper_matrix_index = parent_id;
          }
          BeagleOperation op = {node_id + int_node_count,
                                BEAGLE_OP_NONE,
                                BEAGLE_OP_NONE,
                                upper_partial_index,
                                upper_matrix_index,
                                sister_id,
                                sister_id};
          // Scalers are indexed differently for the upper conditional
          // likelihood. They start at the number of internal nodes + 1 because
          // of the lower conditional likelihoods. Also, in this case the leaves
          // have scalers.
          if (rescaling) {
            // Scaling factors are recomputed every time so we don't read them
            // using destinationScaleRead.
            op.destinationScaleWrite = node_id + 1 + internal_count;
          }
          operations.push_back(op);
        }
      });

  beagleUpdateTransitionMatrices(
      beagle_instance_,
      0,  // eigenIndex
      node_indices.data(),
      gradient_indices.data(),  // firstDerivativeIndices
      nullptr,                  // secondDervativeIndices
      tree.BranchLengths().data(), int_node_count - 1);

  beagleUpdatePartials(beagle_instance_,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       BEAGLE_OP_NONE);  // cumulative scale index

  std::vector<int> category_weight_index = {0};
  std::vector<int> state_frequency_index = {0};
  std::vector<int> cumulative_scale_index = {rescaling_ ? 0 : BEAGLE_OP_NONE};
  std::vector<int> upper_partials_index = {0};
  std::vector<int> node_partial_indices = {0};
  std::vector<int> node_mat_indices = {0};
  std::vector<int> node_deriv_index = {0};
  std::vector<double> gradient(node_count, 0);
  double log_like = 0;

  SizeVectorVector indices_above = tree.Topology()->IdsAbove();
  for (auto &indices : indices_above) {
    // Reverse vector so we have indices from node to root.
    std::reverse(indices.begin(), indices.end());
    // Remove the root scalers.
    indices.pop_back();
  }

  // Actually compute gradient.
  // TODO pass everything by reference
  tree.Topology()->TripleIdPreOrderBifurcating([&](int node_id, int sister_id,
                                                   int) {
    if (node_id != fixed_node_id) {
      double dlogLp;
      upper_partials_index[0] = node_id + int_node_count;
      node_partial_indices[0] = node_id;
      node_mat_indices[0] = node_id;
      node_deriv_index[0] = node_id + int_node_count;

      if (node_id == root_child_id) {
        upper_partials_index[0] = sister_id;
      }
      // parent partial Buffers cannot be a taxon in
      // beagleCalculateEdgeLogLikelihoods
      if (node_partial_indices[0] > upper_partials_index[0]) {
        std::swap(node_partial_indices, upper_partials_index);
      }

      if (rescaling_) {
        beagleResetScaleFactors(beagle_instance_, cumulative_scale_index[0]);
        std::vector<int> scaler_indices(
            static_cast<size_t>(internal_count - 1));
        std::iota(scaler_indices.begin(), scaler_indices.end(), 1);
        // Replace lower scaler index with upper scaler index for nodes between
        // node_index and root.
        int child = node_id;
        for (size_t upper : indices_above[static_cast<size_t>(node_id)]) {
          int int_upper = static_cast<int>(upper);
          int scaler_indices_index = int_upper - int_taxon_count;
          Assert(scaler_indices_index >= 0,
                 "int_upper must be >= taxon count.");
          scaler_indices[static_cast<size_t>(scaler_indices_index)] =
              child + internal_count + 1;
          child = int_upper;
        }
        beagleAccumulateScaleFactors(beagle_instance_, scaler_indices.data(),
                                     static_cast<int>(scaler_indices.size()),
                                     cumulative_scale_index[0]);
      }

      int mysterious_count = 1;  // Not sure what this variable does.
      beagleCalculateEdgeLogLikelihoods(
          beagle_instance_,               // instance number
          upper_partials_index.data(),    // indices of parent partialsBuffers
          node_partial_indices.data(),    // indices of child partialsBuffers
          node_mat_indices.data(),        // transition probability matrices
          node_deriv_index.data(),        // first derivative matrices
          nullptr,                        // second derivative matrices
          category_weight_index.data(),   // pattern weights
          state_frequency_index.data(),   // state frequencies
          cumulative_scale_index.data(),  // scale Factors
          mysterious_count,               // Number of partialsBuffer
          &log_like,                      // destination for log likelihood
          &dlogLp,                        // destination for first derivative
          nullptr);                       // destination for second derivative
      gradient[static_cast<size_t>(node_id)] = dlogLp;
    }
  });

  return {log_like, gradient};
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
