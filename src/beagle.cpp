// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "beagle.hpp"
#include <iostream>
#include <numeric>
#include <string>
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
      std::cerr << "Symbol '" << str[i] << "' not known.\n";
      abort();
    }
  }
  return v;
}

int CreateInstance(int tip_count, int alignment_length,
                   BeagleInstanceDetails *return_info) {
  // Number of partial buffers to create (input) -- internal node count
  // tip_count - 1 for lower partials (internal nodes only)
  // 2*tip_count - 2 for upper partials (every node except the root)
  int partials_buffer_count = 3 * tip_count - 3;
  // Number of compact state representation buffers to create -- for use with
  // setTipStates (input) */
  int compact_buffer_count = tip_count;
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
  // Number of scaling buffers -- can be zero if scaling is not needed
  int scale_buffer_count = 0;
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
  int requirement_flags = 0;

  return beagleCreateInstance(
      tip_count, partials_buffer_count, compact_buffer_count, state_count,
      pattern_count, eigen_buffer_count, matrix_buffer_count, category_count,
      scale_buffer_count, allowed_resources, resource_count, preference_flags,
      requirement_flags, return_info);
}

BeagleInstance CreateInstance(const Alignment &alignment) {
  BeagleInstanceDetails return_info;
  return CreateInstance(static_cast<int>(alignment.SequenceCount()),
                        static_cast<int>(alignment.Length()), &return_info);
}

void SetTipStates(int beagle_instance, const TagStringMap &tag_taxon_map,
                  const Alignment &alignment, const CharIntMap &symbol_table) {
  for (const auto &iter : tag_taxon_map) {
    int taxon_number = static_cast<int>(UnpackFirstInt(iter.first));
    SymbolVector symbols =
        SymbolVectorOf(alignment.at(iter.second), symbol_table);
    beagleSetTipStates(beagle_instance, taxon_number, symbols.data());
  }
}

void PrepareBeagleInstance(
    const BeagleInstance beagle_instance,
    const TreeCollection::TreeCollectionPtr &tree_collection,
    const Alignment &alignment, const CharIntMap &symbol_table) {
  if (tree_collection->TaxonCount() != alignment.SequenceCount()) {
    std::cerr << "The number of tree tips doesn't match the alignment "
                 "sequence count!\n";
    abort();
  }
  SetTipStates(beagle_instance, tree_collection->TagTaxonMap(), alignment,
               symbol_table);
  std::vector<double> pattern_weights(alignment.Length(), 1.);
  beagleSetPatternWeights(beagle_instance, pattern_weights.data());
  // Use uniform rates and weights.
  const double weights[1] = {1.0};
  const double rates[1] = {1.0};
  beagleSetCategoryWeights(beagle_instance, 0, weights);
  beagleSetCategoryRates(beagle_instance, rates);
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

Tree::TreePtr PrepareTreeForLikelihood(Tree::TreePtr tree) {
  if (tree->Children().size() == 3) {
    return tree->Detrifurcate();
  }  // else
  if (tree->Children().size() == 2) {
    return tree;
  }
  // else
  std::cerr << "Tree likelihood calculations should be done on a tree with a "
               "bifurcation or a trifurcation at the root.";
  abort();
}

double LogLikelihood(BeagleInstance beagle_instance, Tree::TreePtr tree) {
  tree = PrepareTreeForLikelihood(tree);
  std::vector<BeagleOperation> operations;
  tree->Topology()->PostOrder([&operations](const Node *node) {
    if (!node->IsLeaf()) {
      assert(node->Children().size() == 2);
      int dest = static_cast<int>(node->Index());
      int child0_index = static_cast<int>(node->Children()[0]->Index());
      int child1_index = static_cast<int>(node->Children()[1]->Index());
      BeagleOperation op = {
          dest,                            // dest
          BEAGLE_OP_NONE, BEAGLE_OP_NONE,  // src and dest scaling
          child0_index,   child0_index,    // src1 and matrix1
          child1_index,   child1_index     // src2 and matrix2
      };
      operations.push_back(op);
    }
  });
  auto branch_count = tree->BranchLengths().size();
  std::vector<int> node_indices(branch_count);
  std::iota(node_indices.begin(), node_indices.end(), 0);
  beagleUpdateTransitionMatrices(beagle_instance,
                                 0,  // eigenIndex
                                 node_indices.data(),
                                 NULL,  // firstDerivativeIndices
                                 NULL,  // secondDervativeIndices
                                 tree->BranchLengths().data(),
                                 static_cast<int>(branch_count));
  beagleUpdatePartials(beagle_instance,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       BEAGLE_OP_NONE);  // cumulative scale index
  double log_like = 0;
  std::vector<int> root_index = {node_indices.back()};
  std::vector<int> category_weight_index = {0};
  std::vector<int> state_frequency_index = {0};
  std::vector<int> cumulative_scale_index = {BEAGLE_OP_NONE};
  int count = 1;
  beagleCalculateRootLogLikelihoods(
      beagle_instance, root_index.data(), category_weight_index.data(),
      state_frequency_index.data(), cumulative_scale_index.data(), count,
      &log_like);
  return log_like;
}

std::vector<double> LogLikelihoods(
    std::vector<BeagleInstance> beagle_instances,
    TreeCollection::TreeCollectionPtr tree_collection) {
  return Parallelize<double>(LogLikelihood, beagle_instances, tree_collection);
}

// Compute first derivative of the log likelihood with respect to each branch
// length, as a vector of first derivatives indexed by node index.
std::vector<double> BranchGradient(BeagleInstance beagle_instance,
                                   Tree::TreePtr tree) {
  tree = PrepareTreeForLikelihood(tree);
  tree->SlideRootPosition();

  std::vector<BeagleOperation> operations;

  size_t branch_count = tree->BranchLengths().size();
  std::vector<int> node_indices(branch_count);
  std::iota(node_indices.begin(), node_indices.end(), 0);

  std::vector<int> gradient_indices(branch_count);
  std::iota(gradient_indices.begin(), gradient_indices.end(),
            tree->BranchLengths().size());

  int int_branch_count = static_cast<int>(branch_count);
  int fixed_node_index =
      static_cast<int>(tree->Topology()->Children()[1]->Index());
  int root_child_index =
      static_cast<int>(tree->Topology()->Children()[0]->Index());

  // Calculate lower partials
  tree->Topology()->BinaryIndexPostOrder(
      [&operations](int node_index, int child0_index, int child1_index) {
        BeagleOperation op = {
            node_index,                      // dest
            BEAGLE_OP_NONE, BEAGLE_OP_NONE,  // src and dest scaling
            child0_index,   child0_index,    // src1 and matrix1
            child1_index,   child1_index     // src2 and matrix2
        };
        operations.push_back(op);
      });

  // Calculate upper partials
  tree->Topology()->TripleIndexPreOrderBifurcating(
      [&operations, &root_child_index, &int_branch_count, &fixed_node_index](
          int parent_index, int sister_index, int node_index) {
        if (node_index != root_child_index && node_index != fixed_node_index) {
          int upper_partial_index;
          int upper_matrix_index;
          if (parent_index == root_child_index) {
            upper_matrix_index = static_cast<int>(root_child_index);
            upper_partial_index = static_cast<int>(fixed_node_index);
          } else if (parent_index == fixed_node_index) {
            upper_matrix_index = static_cast<int>(root_child_index);
            upper_partial_index = static_cast<int>(root_child_index);
          } else {
            upper_partial_index = parent_index + int_branch_count;
            upper_matrix_index = parent_index;
          }
          BeagleOperation op = {node_index + int_branch_count,
                                BEAGLE_OP_NONE,
                                BEAGLE_OP_NONE,
                                upper_partial_index,
                                upper_matrix_index,
                                sister_index,
                                sister_index};
          operations.push_back(op);
        }
      });

  beagleUpdateTransitionMatrices(
      beagle_instance,
      0,  // eigenIndex
      node_indices.data(),
      gradient_indices.data(),  // firstDerivativeIndices
      NULL,                     // secondDervativeIndices
      tree->BranchLengths().data(), int_branch_count);
  beagleUpdatePartials(beagle_instance,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       BEAGLE_OP_NONE);  // cumulative scale index

  std::vector<int> category_weight_index = {0};
  std::vector<int> state_frequency_index = {0};
  std::vector<int> cumulative_scale_index = {BEAGLE_OP_NONE};
  int count = 1;
  std::vector<int> upper_partials_index = {0};
  std::vector<int> node_partial_indices = {0};
  std::vector<int> node_mat_indices = {0};
  std::vector<int> node_deriv_index = {0};
  std::vector<double> derivatives(branch_count, 0);

  // Actually compute gradient.
  tree->Topology()->TripleIndexPreOrderBifurcating(
      [&](int, int sister_index, int node_index) {
        if (node_index != fixed_node_index) {
          double dlogLp;
          double log_like = 0;
          upper_partials_index[0] = node_index + int_branch_count;
          node_partial_indices[0] = node_index;
          node_mat_indices[0] = node_index;
          node_deriv_index[0] = node_index + int_branch_count;
          if (node_index == root_child_index) {
            upper_partials_index[0] = sister_index;
          }
          // parent partial Buffers cannot be a taxon in
          // beagleCalculateEdgeLogLikelihoods
          if (node_partial_indices[0] > upper_partials_index[0]) {
            int temp = node_partial_indices[0];
            node_partial_indices[0] = upper_partials_index[0];
            upper_partials_index[0] = temp;
          }
          beagleCalculateEdgeLogLikelihoods(
              beagle_instance,              // instance number
              upper_partials_index.data(),  // indices of parent partialsBuffers
              node_partial_indices.data(),  // indices of child partialsBuffers
              node_mat_indices.data(),      // transition probability matrices
              node_deriv_index.data(),      // first derivative matrices
              NULL,                         // second derivative matrices
              category_weight_index.data(),   // pattern weights
              state_frequency_index.data(),   // state frequencies
              cumulative_scale_index.data(),  // scale Factors
              count,                          // Number of partialsBuffer
              &log_like,                      // destination for log likelihood
              &dlogLp,  // destination for first derivative  # derivative code
              NULL);    // destination for second derivative
          derivatives[static_cast<size_t>(node_index)] = dlogLp;
        }
      });

  return derivatives;
}

std::vector<std::vector<double>> BranchGradients(
    std::vector<BeagleInstance> beagle_instances,
    TreeCollection::TreeCollectionPtr tree_collection) {
  return Parallelize<std::vector<double>>(BranchGradient, beagle_instances,
                                          tree_collection);
}

}  // namespace beagle
