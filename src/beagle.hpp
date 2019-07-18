// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BEAGLE_HPP_
#define SRC_BEAGLE_HPP_

#include <string>
#include <vector>
#include "alignment.hpp"
#include "intpack.hpp"
#include "libhmsbeagle/beagle.h"
#include "tree_collection.hpp"
#include "typedefs.hpp"

namespace beagle {

typedef int BeagleInstance;

CharIntMap GetSymbolTable();
SymbolVector SymbolVectorOf(const std::string &str,
                            const CharIntMap &symbol_table);

int CreateInstance(int tip_count, int alignment_length,
                   BeagleInstanceDetails *return_info);
BeagleInstance CreateInstance(const Alignment &alignment);

void SetTipStates(int beagle_instance, const TagStringMap &tag_taxon_map,
                  const Alignment &alignment, const CharIntMap &symbol_table);
void PrepareBeagleInstance(
    const BeagleInstance beagle_instance,
    const TreeCollection::TreeCollectionPtr &tree_collection,
    const Alignment &alignment, const CharIntMap &symbol_table);
void SetJCModel(BeagleInstance beagle_instance);

double BifurcatingTreeLogLikelihood(Tree::TreePtr tree,
                                    BeagleInstance beagle_instance);
double TreeLogLikelihood(Tree::TreePtr tree, BeagleInstance beagle_instance);
std::vector<double> TreeLogLikelihoods(
    BeagleInstance beagle_instance,
    TreeCollection::TreeCollectionPtr tree_collection);
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
