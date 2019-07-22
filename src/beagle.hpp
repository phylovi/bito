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
std::vector<T> Parallelize(
    std::function<T(BeagleInstance, Tree::TreePtr)> f,
    std::vector<BeagleInstance> beagle_instances,
    TreeCollection::TreeCollectionPtr tree_collection) {
  if (beagle_instances.size() == 0) {
    std::cerr << "Please add some BEAGLE instances that can be used for "
                 "computation.\n";
    abort();
  }
  std::vector<T> results(tree_collection->TreeCount());
  std::queue<BeagleInstance> instance_queue;
  for (auto instance : beagle_instances) {
    instance_queue.push(instance);
  }
  std::queue<size_t> tree_number_queue;
  for (size_t i = 0; i < tree_collection->TreeCount(); i++) {
    tree_number_queue.push(i);
  }
  TaskProcessor<BeagleInstance, size_t> task_processor(
      instance_queue, tree_number_queue,
      [&results, &tree_collection, &f](BeagleInstance beagle_instance,
                                       size_t tree_number) {
        results[tree_number] =
            f(beagle_instance, tree_collection->GetTree(tree_number));
      });
  return results;
}

double LogLikelihood(BeagleInstance beagle_instance, Tree::TreePtr tree);
std::vector<double> LogLikelihoods(
    std::vector<BeagleInstance> beagle_instances,
    TreeCollection::TreeCollectionPtr tree_collection);

std::vector<double> BranchGradient(BeagleInstance beagle_instance,
                                   Tree::TreePtr tree);
std::vector<std::vector<double>> BranchGradients(
    std::vector<BeagleInstance> beagle_instances,
    TreeCollection::TreeCollectionPtr tree_collection);

}  // namespace beagle

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Beagle") {
  CharIntMap symbol_table = beagle::GetSymbolTable();
  SymbolVector symbol_vector =
      beagle::SymbolVectorOf("-tgcaTGCA", symbol_table);
  SymbolVector correct_symbol_vector = {4, 3, 2, 1, 0, 3, 2, 1, 0};
  CHECK_EQ(symbol_vector, correct_symbol_vector);

  // The real tests are in libsbn.hpp, where we have access to tree parsing.
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_BEAGLE_HPP_
