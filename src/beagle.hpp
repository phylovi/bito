// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BEAGLE_HPP_
#define SRC_BEAGLE_HPP_

#include <numeric>
#include <string>
#include "alignment.hpp"
#include "intpack.hpp"
#include "libhmsbeagle/beagle.h"
#include "tree_collection.hpp"
#include "typedefs.hpp"

namespace beagle {

typedef int BeagleInstance;

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
  int partials_buffer_count = tip_count - 1;
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
  // Number of transition matrix buffers (input) -- one per edge
  int matrix_buffer_count = 2 * tip_count - 1;
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
  long preference_flags = 0;
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
  BeagleInstanceDetails *return_info = new BeagleInstanceDetails();
  // Not worrying about freeing this return_info.
  return CreateInstance(static_cast<int>(alignment.SequenceCount()),
                        static_cast<int>(alignment.Length()), return_info);
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

double BifurcatingTreeLogLikelihood(Tree::TreePtr tree,
                                    BeagleInstance beagle_instance) {
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

double TreeLogLikelihood(Tree::TreePtr tree, BeagleInstance beagle_instance) {
  if (tree->Children().size() == 2) {
    return BifurcatingTreeLogLikelihood(tree, beagle_instance);
  }
  // else
  if (tree->Children().size() == 3) {
    return BifurcatingTreeLogLikelihood(tree->Detrifurcate(), beagle_instance);
  }
  // else
  std::cerr << "Tree likelihood calculations should be done on a tree with a "
               "bifurcation or a trifurcation at the root.";
  abort();
}

std::vector<double> TreeLogLikelihoods(
    BeagleInstance beagle_instance,
    TreeCollection::TreeCollectionPtr tree_collection) {
  std::vector<double> results;
  for (const auto &tree : tree_collection->Trees()) {
    results.push_back(TreeLogLikelihood(tree, beagle_instance));
  }
  return results;
}

}  // namespace beagle

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Beagle") {
  CharIntMap symbol_table = beagle::GetSymbolTable();
  SymbolVector symbol_vector =
      beagle::SymbolVectorOf("-tgcaTGCA", symbol_table);
  SymbolVector correct_symbol_vector = {4, 3, 2, 1, 0, 3, 2, 1, 0};
  CHECK_EQ(symbol_vector, correct_symbol_vector);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_BEAGLE_HPP_
