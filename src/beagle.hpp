// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BEAGLE_HPP_
#define SRC_BEAGLE_HPP_

#include <functional>
#include <queue>
#include <string>
#include <utility>
#include <vector>
#include "alignment.hpp"
#include "intpack.hpp"
#include "libhmsbeagle/beagle.h"
#include "site_pattern.hpp"
#include "sugar.hpp"
#include "task_processor.hpp"
#include "tree_collection.hpp"

namespace beagle {

typedef int BeagleInstance;

int CreateInstance(int tip_count, int alignment_length,
                   BeagleInstanceDetails *return_info);
BeagleInstance CreateInstance(const SitePattern &site_pattern);

void PrepareBeagleInstance(const BeagleInstance beagle_instance,
                           const TreeCollection &tree_collection,
                           const SitePattern &site_pattern);

void SetJCModel(BeagleInstance beagle_instance);
void SetSubstModel(BeagleInstance beagle_instance,
                   std::vector<double> freqs,
                   std::vector<double> evec,
                   std::vector<double> ivec,
                   std::vector<double> eval);


template <typename T>
std::vector<T> Parallelize(
    std::function<T(BeagleInstance, const Tree &, bool rescaling)> f,
    std::vector<BeagleInstance> beagle_instances,
    const TreeCollection &tree_collection, bool rescaling) {
  if (beagle_instances.size() == 0) {
    Failwith(
        "Please add some BEAGLE instances that can be used for computation.");
  }
  std::vector<T> results(tree_collection.TreeCount());
  std::queue<BeagleInstance> instance_queue;
  for (auto instance : beagle_instances) {
    instance_queue.push(instance);
  }
  std::queue<size_t> tree_number_queue;
  for (size_t i = 0; i < tree_collection.TreeCount(); i++) {
    tree_number_queue.push(i);
  }
  TaskProcessor<BeagleInstance, size_t> task_processor(
      std::move(instance_queue), std::move(tree_number_queue),
      [&results, &tree_collection, &f, &rescaling](
          BeagleInstance beagle_instance, size_t tree_number) {
        results[tree_number] =
            f(beagle_instance, tree_collection.GetTree(tree_number), rescaling);
      });
  return results;
}

double LogLikelihood(BeagleInstance beagle_instance, const Tree &tree,
                     bool rescaling);
std::vector<double> LogLikelihoods(std::vector<BeagleInstance> beagle_instances,
                                   const TreeCollection &tree_collection,
                                   bool rescaling);

std::pair<double, std::vector<double>> BranchGradient(
    BeagleInstance beagle_instance, const Tree &tree, bool rescaling);
std::vector<std::pair<double, std::vector<double>>> BranchGradients(
    std::vector<BeagleInstance> beagle_instances,
    const TreeCollection &tree_collection, bool rescaling);

std::vector<double> Gradient(BeagleInstance beagle_instance, const Tree &tree, bool rescaling);

double ParameterGradient(BeagleInstance beagle_instance, const Tree &tree, bool rescaling, std::vector<double> q_differential,
                   std::vector<double> freqs,
                   std::vector<double> evec,
                   std::vector<double> ivec,
                   std::vector<double> eval);
}  // namespace beagle

// The tests are in libsbn.hpp, where we have access to tree parsing.

#endif  // SRC_BEAGLE_HPP_
