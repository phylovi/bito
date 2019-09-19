// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "beagle.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace beagle {

// DNA assumption here.
CharIntMap GetSymbolTable() {
  CharIntMap table({{'A', 0},
                    {'C', 1},
                    {'G', 2},
                    {'T', 3},
                    {'a', 0},
                    {'c', 1},
                    {'g', 2},
                    {'t', 3},
                    {'-', 4}});
  return table;
}

SymbolVector SymbolVectorOf(const std::string &str,
                            const CharIntMap &symbol_table) {
  SymbolVector v(str.size());
  for (size_t i = 0; i < str.size(); i++) {
    auto search = symbol_table.find(str[i]);
    if (search != symbol_table.end()) {
      v[i] = search->second;
    } else {
      Failwith("Symbol '" + std::to_string(str[i]) + "' not known.\n");
    }
  }
  return v;
}

int CreateInstance(int tip_count, int alignment_length,
                   BeagleInstanceDetails *return_info) {
  bool use_tip_states = false;
  // Number of partial buffers to create (input):
  // tip_count - 1 for lower partials (internal nodes only)
  // 2*tip_count - 1 for upper partials (every node)
  int partials_buffer_count = 3 * tip_count - 2;
  if (!use_tip_states) partials_buffer_count += tip_count;
  // Number of compact state representation buffers to create -- for use with
  // setTipStates (input) */
  int compact_buffer_count = (use_tip_states ? tip_count : 0);
  // DNA assumption here.
  int state_count = 4;
  // Number of site patterns to be handled by the instance (input) -- not
  // compressed in this case
  int pattern_count = alignment_length;
  // Number of eigen-decomposition buffers to allocate (input)
  int eigen_buffer_count = 1;
  // Number of transition matrix buffers (input) -- two per edge
  int matrix_buffer_count = 2 * (2 * tip_count - 1);
  // Number of rate categories
  int category_count = 1;
  // Number of scaling buffers -- 1 buffer per partial buffer and 1 more
  // for accumulating scale factors in position 0.
  int scale_buffer_count = partials_buffer_count + 1;
  // List of potential resources on which this instance is allowed (input,
  // NULL implies no restriction
  int *allowed_resources = nullptr;
  // Length of resourceList list (input) -- not needed to use the default
  // hardware config
  int resource_count = 0;
  // Bit-flags indicating preferred implementation charactertistics, see
  // BeagleFlags (input)
  int64_t preference_flags = 0;
  // Bit-flags indicating required implementation characteristics, see
  // BeagleFlags (input)
  int requirement_flags = BEAGLE_FLAG_SCALING_MANUAL;

  return beagleCreateInstance(
      tip_count, partials_buffer_count, compact_buffer_count, state_count,
      pattern_count, eigen_buffer_count, matrix_buffer_count, category_count,
      scale_buffer_count, allowed_resources, resource_count, preference_flags,
      requirement_flags, return_info);
}

BeagleInstance CreateInstance(const SitePattern &site_pattern) {
  BeagleInstanceDetails return_info;
  return CreateInstance(static_cast<int>(site_pattern.SequenceCount()),
                        static_cast<int>(site_pattern.PatternCount()),
                        &return_info);
}

void PrepareBeagleInstance(const BeagleInstance beagle_instance,
                           const TreeCollection &tree_collection,
                           const SitePattern &site_pattern) {
  if (tree_collection.TaxonCount() != site_pattern.SequenceCount()) {
    Failwith(
        "The number of tree tips doesn't match the alignment sequence count.");
  }
  // Use uniform rates and weights.
  const double weights[1] = {1.0};
  const double rates[1] = {1.0};
  const std::vector<SymbolVector> patterns = site_pattern.GetPatterns();
  bool use_tip_states = false;
  if (use_tip_states) {
    int taxon_number = 0;
    for (const auto &pattern : site_pattern.GetPatterns()) {
      beagleSetTipStates(beagle_instance, taxon_number++, pattern.data());
    }
  } else {
    // In case we want to use tip partials instead of tip states
    for (int i = 0; i < site_pattern.GetPatterns().size(); i++) {
      beagleSetTipPartials(beagle_instance, i,
                           site_pattern.GetPartials(i).data());
    }
  }
  beagleSetPatternWeights(beagle_instance, site_pattern.GetWeights().data());
  beagleSetCategoryWeights(beagle_instance, 0, weights);
  beagleSetCategoryRates(beagle_instance, rates);

  // Set the preorder partials of root node (state frequencies)
  // For JC69 the frequencies are fixed but for other models the code below will
  // have to move where preorder partials are calculated
  std::vector<double> state_frequencies(site_pattern.PatternCount() * 4, 0.25);
  int int_node_count = static_cast<int>(2 * site_pattern.SequenceCount() - 1);
  int root_id = static_cast<int>(tree_collection.GetTree(0).Id()) +
                1;  // + because trees are detrifurcated
  beagleSetPartials(beagle_instance, root_id + int_node_count,
                    state_frequencies.data());
}

void SetJCModel(BeagleInstance beagle_instance) {
  std::vector<double> freqs(4, 0.25);
  beagleSetStateFrequencies(beagle_instance, 0, freqs.data());
  // an eigen decomposition for the JC69 model
  std::vector<double> evec = {1.0, 2.0, 0.0, 0.5,  1.0, -2.0, 0.5,  0.0,
                              1.0, 2.0, 0.0, -0.5, 1.0, -2.0, -0.5, 0.0};

  std::vector<double> ivec = {0.25,  0.25,   0.25, 0.25, 0.125, -0.125,
                              0.125, -0.125, 0.0,  1.0,  0.0,   -1.0,
                              1.0,   0.0,    -1.0, 0.0};

  std::vector<double> eval = {0.0, -1.3333333333333333, -1.3333333333333333,
                              -1.3333333333333333};
  beagleSetEigenDecomposition(beagle_instance, 0, evec.data(), ivec.data(),
                              eval.data());
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

double LogLikelihood(BeagleInstance beagle_instance, const Tree &in_tree,
                     bool rescaling) {
  beagleResetScaleFactors(beagle_instance, 0);
  auto tree = PrepareTreeForLikelihood(in_tree);
  auto node_count = tree.BranchLengths().size();
  int int_taxon_count = static_cast<int>(tree.LeafCount());
  std::vector<BeagleOperation> operations;
  tree.Topology()->PostOrder(
      [&operations, int_taxon_count, rescaling](const Node *node) {
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
  beagleUpdateTransitionMatrices(beagle_instance,
                                 0,  // eigenIndex
                                 node_indices.data(),
                                 NULL,  // firstDerivativeIndices
                                 NULL,  // secondDervativeIndices
                                 tree.BranchLengths().data(),
                                 static_cast<int>(node_count - 1));

  // This is the entry of scaleBuffer in which we store accumulated factors.
  int cumululative_index = rescaling ? 0 : BEAGLE_OP_NONE;
  beagleUpdatePartials(beagle_instance,
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
      beagle_instance, root_id.data(), category_weight_index.data(),
      state_frequency_index.data(), cumulative_scale_index.data(),
      mysterious_count, &log_like);
  return log_like;
}

std::vector<double> LogLikelihoods(std::vector<BeagleInstance> beagle_instances,
                                   const TreeCollection &tree_collection,
                                   bool rescaling) {
  return Parallelize<double>(LogLikelihood, beagle_instances, tree_collection,
                             rescaling);
}

// Compute first derivative of the log likelihood with respect to each branch
// length, as a vector of first derivatives indexed by node id.
std::pair<double, std::vector<double>> BranchGradient(
    BeagleInstance beagle_instance, const Tree &in_tree, bool rescaling) {
  beagleResetScaleFactors(beagle_instance, 0);
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
      [&operations, int_taxon_count, rescaling](int node_id, int child0_id,
                                                int child1_id) {
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
      [&operations, &root_child_id, &fixed_node_id, rescaling, internal_count,
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
      beagle_instance,
      0,  // eigenIndex
      node_indices.data(),
      gradient_indices.data(),  // firstDerivativeIndices
      NULL,                     // secondDervativeIndices
      tree.BranchLengths().data(), int_node_count - 1);

  beagleUpdatePartials(beagle_instance,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       BEAGLE_OP_NONE);  // cumulative scale index

  std::vector<int> category_weight_index = {0};
  std::vector<int> state_frequency_index = {0};
  std::vector<int> cumulative_scale_index = {rescaling ? 0 : BEAGLE_OP_NONE};
  int mysterious_count = 1;
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

      if (rescaling) {
        beagleResetScaleFactors(beagle_instance, cumulative_scale_index[0]);
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
        beagleAccumulateScaleFactors(beagle_instance, scaler_indices.data(),
                                     static_cast<int>(scaler_indices.size()),
                                     cumulative_scale_index[0]);
      }

      beagleCalculateEdgeLogLikelihoods(
          beagle_instance,                // instance number
          upper_partials_index.data(),    // indices of parent partialsBuffers
          node_partial_indices.data(),    // indices of child partialsBuffers
          node_mat_indices.data(),        // transition probability matrices
          node_deriv_index.data(),        // first derivative matrices
          NULL,                           // second derivative matrices
          category_weight_index.data(),   // pattern weights
          state_frequency_index.data(),   // state frequencies
          cumulative_scale_index.data(),  // scale Factors
          mysterious_count,               // Number of partialsBuffer
          &log_like,                      // destination for log likelihood
          &dlogLp,                        // destination for first derivative
          NULL);                          // destination for second derivative
      gradient[static_cast<size_t>(node_id)] = dlogLp;
    }
  });

  return {log_like, gradient};
}

std::vector<std::pair<double, std::vector<double>>> BranchGradients(
    std::vector<BeagleInstance> beagle_instances,

    const TreeCollection &tree_collection, bool rescaling) {
  return Parallelize<std::pair<double, std::vector<double>>>(
      BranchGradient, beagle_instances, tree_collection, rescaling);
}

std::vector<std::pair<double, std::vector<double>>> NewBranchGradients(
    std::vector<BeagleInstance> beagle_instances,

    const TreeCollection &tree_collection, bool rescaling) {
  return Parallelize<std::pair<double, std::vector<double>>>(
      NewBranchGradient, beagle_instances, tree_collection, rescaling);
}

/// Compute first derivative of the log likelihood with respect to each branch
// length, as a vector of first derivatives indexed by node id.
std::pair<double, std::vector<double>> NewBranchGradient(
    BeagleInstance beagle_instance, const Tree &in_tree, bool rescaling) {
  beagleResetScaleFactors(beagle_instance, 0);
  auto tree = PrepareTreeForLikelihood(in_tree);
  tree.SlideRootPosition();

  size_t node_count = tree.BranchLengths().size();

  int int_node_count = static_cast<int>(node_count);
  int int_taxon_count = static_cast<int>(tree.LeafCount());
  int internal_count = int_taxon_count - 1;
  std::vector<int> node_indices(node_count - 1);
  std::iota(node_indices.begin(), node_indices.end(), 0);
  std::vector<int> derivative_matrix_indices(int_node_count - 1);
  std::iota(derivative_matrix_indices.begin(), derivative_matrix_indices.end(),
            node_count - 1);
  std::vector<BeagleOperation> operations;
  std::vector<BeagleOperation> pre_operations;

  int root_id = static_cast<int>(tree.Topology()->Id());

  // Calculate post-order partials
  tree.Topology()->BinaryIdPostOrder([&operations, int_taxon_count, rescaling](
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

  beagleUpdateTransitionMatrices(
      beagle_instance,
      0,  // eigenIndex
      node_indices.data(),
      derivative_matrix_indices.data(),  // firstDerivativeIndices
      NULL,                              // secondDervativeIndices
      tree.BranchLengths().data(), int_node_count - 1);

  beagleUpdatePartials(beagle_instance, operations.data(),
                       static_cast<int>(operations.size()),
                       BEAGLE_OP_NONE);  // cumulative scale index

  // Calculate pre-order partials
  // TODO: Set root pre-order partials equal to root state frequencies

  tree.Topology()->TripleIdPreOrderBifurcating([&pre_operations, &root_id,
                                                &int_node_count, &rescaling,
                                                &internal_count](
      int node_id, int sister_id, int parent_id) {
    if (node_id != root_id) {
      BeagleOperation op = {
          node_id + int_node_count,  // dest pre-order partial of current node
          BEAGLE_OP_NONE,
          BEAGLE_OP_NONE,
          parent_id + int_node_count,  // pre-order partial parent
          node_id,                     // matrices of current node
          sister_id,                   // post-order partial of sibling
          sister_id                    // matrices of sibling
      };
      // Scalers are indexed differently for the upper conditional
      // likelihood. They start at the number of internal nodes + 1 because
      // of the lower conditional likelihoods. Also, in this case the leaves
      // have scalers.
      if (rescaling) {
        // Scaling factors are recomputed every time so we don't read them
        // using destinationScaleRead.
        op.destinationScaleWrite = node_id + 1 + internal_count;
      }
      pre_operations.push_back(op);
    }
  });

  beagleUpdatePrePartials(beagle_instance, pre_operations.data(),
                          static_cast<int>(pre_operations.size()),
                          BEAGLE_OP_NONE);  // cumulative scale index

  std::vector<double> gradient(node_count, 0);
  std::vector<int> category_weight_index = {0};
  std::vector<int> post_buffer_indices(int_node_count - 1);
  std::vector<int> pre_buffer_indices(int_node_count - 1);
  std::iota(post_buffer_indices.begin(), post_buffer_indices.end(), 0);
  std::iota(pre_buffer_indices.begin(), pre_buffer_indices.end(),
            int_node_count);

  beagleCalculateEdgeDerivatives(
      beagle_instance,
      post_buffer_indices.data(),        // list of post order buffer indices
      pre_buffer_indices.data(),         // list of pre order buffer indices
      derivative_matrix_indices.data(),  // differential Q matrix indices
      category_weight_index.data(),      // category weights indices
      int_node_count - 1,                // number of edges
      NULL,                              // derivative-per-site output array
      gradient.data(),  // sum of derivatives across sites output array
      NULL);            // sum of squared derivatives output array

  return {0, gradient};
}

}  // namespace beagle
