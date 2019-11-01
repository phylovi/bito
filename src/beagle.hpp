// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_BEAGLE_HPP_
#define SRC_BEAGLE_HPP_

#include <functional>
#include <memory>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "alignment.hpp"
#include "clock_model.hpp"
#include "intpack.hpp"
#include "libhmsbeagle/beagle.h"
#include "site_model.hpp"
#include "site_pattern.hpp"
#include "substitution_model.hpp"
#include "sugar.hpp"
#include "task_processor.hpp"
#include "tree_collection.hpp"

using BeagleInstance = int;

// TODO I feel like the presence of rescaling here begs for having a likelihood
// computation method that knows if rescaling is being applied.
template <typename T>
std::vector<T> Parallelize(
    std::function<T(BeagleInstance, const Tree &, bool)> f,
    const std::vector<BeagleInstance> &beagle_instances,
    const TreeCollection &tree_collection, bool rescaling) {
  if (beagle_instances.empty()) {
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

class Engine {
 public:
  Engine(std::unique_ptr<SubstitutionModel> substitution_model,
         std::unique_ptr<SiteModel> site_model,
         std::unique_ptr<ClockModel> clock_model, size_t thread_count,
         const SitePattern &site_pattern);

  void UpdateEigenDecompositionModel();

  void UpdateSiteModel();

  std::vector<double> LogLikelihoods(const TreeCollection &tree_collection);

  std::vector<std::pair<double, std::vector<double>>> BranchGradients(
      const TreeCollection &tree_collection);

 private:
  void CreateInstances(const SitePattern &site_pattern);
  void SetTipStates(const SitePattern &site_pattern);

  std::vector<int> beagle_instances_;
  std::unique_ptr<SubstitutionModel> substitution_model_;
  std::unique_ptr<SiteModel> site_model_;
  std::unique_ptr<ClockModel> clock_model_;
  size_t thread_count_;
  int pattern_count_;
  bool rescaling_;
};

// The tests are in libsbn.hpp, where we have access to tree parsing.

#endif  // SRC_BEAGLE_HPP_
