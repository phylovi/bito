// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BEAGLE_HPP_
#define SRC_BEAGLE_HPP_

#include <functional>
#include <queue>
#include <string>
#include <vector>
#include "alignment.hpp"
#include "intpack.hpp"
#include "libhmsbeagle/beagle.h"
#include "task_processor.hpp"
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

template <typename T>
T Parallelize(std::function<T(Tree::TreePtr, beagle::BeagleInstance)> f,
              TreeCollection::TreeCollectionPtr tree_collection,
              std::vector<beagle::BeagleInstance> beagle_instances) {
  if (beagle_instances.size() == 0) {
    std::cerr << "Please add some BEAGLE instances that can be used for "
                 "computation.\n";
    abort();
  }
  std::vector<double> results(tree_collection->TreeCount());
  std::queue<beagle::BeagleInstance> instance_queue;
  for (auto instance : beagle_instances) {
    instance_queue.push(instance);
  }
  std::queue<size_t> tree_number_queue;
  for (size_t i = 0; i < tree_collection->TreeCount(); i++) {
    tree_number_queue.push(i);
  }
  TaskProcessor<beagle::BeagleInstance, size_t> task_processor(
      instance_queue, tree_number_queue,
      [&results, &tree_collection, &f](beagle::BeagleInstance beagle_instance,
                                       size_t tree_number) {
        results[tree_number] =
            f(tree_collection->GetTree(tree_number), beagle_instance);
      });
  return results;
}

double LogLikelihood(Tree::TreePtr tree, BeagleInstance beagle_instance);

std::vector<double> BranchGradients(Tree::TreePtr tree,
                                    BeagleInstance beagle_instance);

double LogLikelihood(Tree::TreePtr tree, BeagleInstance beagle_instance);
std::vector<double> LogLikelihoods(
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
