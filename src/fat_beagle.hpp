// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_FAT_BEAGLE_HPP_
#define SRC_FAT_BEAGLE_HPP_

#include <memory>
#include <queue>
#include <utility>
#include <vector>
#include "libhmsbeagle/beagle.h"
#include "phylo_model.hpp"
#include "site_pattern.hpp"
#include "task_processor.hpp"
#include "tree_collection.hpp"

class FatBeagle {
 public:
  // This constructor makes the beagle_instance_;
  FatBeagle(const PhyloModelSpecification &specification,
            const SitePattern &site_pattern);
  ~FatBeagle();
  // Delete (copy + move) x (constructor + assignment)
  FatBeagle(const FatBeagle &) = delete;
  FatBeagle(const FatBeagle &&) = delete;
  FatBeagle &operator=(const FatBeagle &) = delete;
  FatBeagle &operator=(const FatBeagle &&) = delete;

  const BlockSpecification &GetBlockSpecification() const;

  void SetParameters(const EigenVectorXdRef param_vector);
  void SetRescaling();

  double LogLikelihood(const Tree &tree) const;
  std::pair<double, std::vector<double>> BranchGradient(const Tree &tree) const;

  // We can pass these static methods to FatBeagleParallelize.
  static double StaticLogLikelihood(FatBeagle *fat_beagle, const Tree &in_tree);
  static std::pair<double, std::vector<double>> StaticBranchGradient(
      FatBeagle *fat_beagle, const Tree &in_tree);

 private:
  using BeagleInstance = int;

  BeagleInstance CreateInstance(const SitePattern &site_pattern);
  void SetTipStates(const SitePattern &site_pattern);
  void UpdateSiteModelInBeagle();
  void UpdateSubstitutionModelInBeagle();
  void UpdatePhyloModelInBeagle();

  std::unique_ptr<PhyloModel> phylo_model_;
  bool rescaling_;
  BeagleInstance beagle_instance_;
  int pattern_count_;
};

template <typename T>
std::vector<T> FatBeagleParallelize(
    std::function<T(FatBeagle *, const Tree &)> f,
    const std::vector<std::unique_ptr<FatBeagle>> &fat_beagles,
    const TreeCollection &tree_collection, EigenMatrixXdRef param_matrix) {
  if (fat_beagles.empty()) {
    Failwith("Please add some FatBeagles that can be used for computation.");
  }
  std::vector<T> results(tree_collection.TreeCount());
  std::queue<FatBeagle *> fat_beagle_queue;
  for (const auto &fat_beagle : fat_beagles) {
    Assert(fat_beagle != nullptr, "Got a fat_beagle nullptr!");
    fat_beagle_queue.push(fat_beagle.get());
  }
  std::queue<size_t> tree_number_queue;
  for (size_t i = 0; i < tree_collection.TreeCount(); i++) {
    tree_number_queue.push(i);
  }
  Assert(tree_collection.TreeCount() == param_matrix.rows(),
         "We param_matrix needs as many rows as we have trees.");
  TaskProcessor<FatBeagle *, size_t> task_processor(
      std::move(fat_beagle_queue), std::move(tree_number_queue),
      [&results, &tree_collection, &param_matrix, &f](FatBeagle *fat_beagle,
                                                      size_t tree_number) {
        fat_beagle->SetParameters(param_matrix.row(tree_number));
        results[tree_number] =
            f(fat_beagle, tree_collection.GetTree(tree_number));
      });
  return results;
}

// Tests live in libsbn.hpp.
#endif  // SRC_FAT_BEAGLE_HPP_
