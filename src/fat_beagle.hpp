// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_FAT_BEAGLE_HPP_
#define SRC_FAT_BEAGLE_HPP_

#include <memory>
#include <queue>
#include <utility>
#include <vector>
#include "beagle_accessories.hpp"
#include "phylo_model.hpp"
#include "rooted_tree_collection.hpp"
#include "site_pattern.hpp"
#include "task_processor.hpp"
#include "tree_collection.hpp"
#include "unrooted_tree_collection.hpp"

class FatBeagle {
 public:
  using PackedBeagleFlags = long;

  // This constructor makes the beagle_instance_;
  FatBeagle(const PhyloModelSpecification &specification,
            const SitePattern &site_pattern,
            const PackedBeagleFlags beagle_preference_flags, bool use_tip_states);
  ~FatBeagle();
  // Delete (copy + move) x (constructor + assignment) because FatBeagle manages an
  // external resource (a BEAGLE instance).
  FatBeagle(const FatBeagle &) = delete;
  FatBeagle(const FatBeagle &&) = delete;
  FatBeagle &operator=(const FatBeagle &) = delete;
  FatBeagle &operator=(const FatBeagle &&) = delete;

  const BlockSpecification &GetPhyloModelBlockSpecification() const;
  const PackedBeagleFlags &GetBeagleFlags() const { return beagle_flags_; };

  void SetParameters(const EigenVectorXdRef param_vector);
  void SetRescaling(const bool rescaling) { rescaling_ = rescaling; }

  double LogLikelihood(const UnrootedTree &tree) const;
  double LogLikelihood(const RootedTree &tree) const;
  // Compute first derivative of the log likelihood with respect to each branch
  // length, as a vector of first derivatives indexed by node id.
  std::pair<double, std::vector<double>> BranchGradient(const UnrootedTree &tree) const;
  std::pair<double, std::vector<double>> BranchGradient(const RootedTree &tree) const;

  // We can pass these static methods to FatBeagleParallelize.
  static double StaticLogLikelihood(FatBeagle *fat_beagle, const UnrootedTree &in_tree);
  static double StaticRootedLogLikelihood(FatBeagle *fat_beagle,
                                          const RootedTree &in_tree);
  static std::pair<double, std::vector<double>> StaticBranchGradient(
      FatBeagle *fat_beagle, const UnrootedTree &in_tree);
  static std::pair<double, std::vector<double>> StaticRootedBranchGradient(
      FatBeagle *fat_beagle, const RootedTree &in_tree);

 private:
  using BeagleInstance = int;
  using BeagleOperationVector = std::vector<BeagleOperation>;

  std::unique_ptr<PhyloModel> phylo_model_;
  bool rescaling_;
  BeagleInstance beagle_instance_;
  PackedBeagleFlags beagle_flags_;
  int pattern_count_;
  bool use_tip_states_;

  std::pair<BeagleInstance, PackedBeagleFlags> CreateInstance(
      const SitePattern &site_pattern, PackedBeagleFlags beagle_preference_flags);
  void SetTipStates(const SitePattern &site_pattern);
  void SetTipPartials(const SitePattern &site_pattern);
  void UpdateSiteModelInBeagle();
  void UpdateSubstitutionModelInBeagle();
  void UpdatePhyloModelInBeagle();

  std::pair<double, std::unordered_map<std::string, std::vector<double>>> Gradient(
      const RootedTree &tree) const;

  std::pair<double, std::unordered_map<std::string, std::vector<double>>> Gradient(
      const UnrootedTree &in_tree) const;

  double LogLikelihoodInternals(const Node::NodePtr topology,
                                const std::vector<double> &branch_lengths) const;
  std::pair<double, std::vector<double>> BranchGradientInternals(
      const Node::NodePtr topology, const std::vector<double> &branch_lengths) const;

  void UpdateBeagleTransitionMatrices(
      const BeagleAccessories &baBranchGradientInternals,
      const std::vector<double> &branch_lengths,
      const int *const gradient_indices_ptr) const;
  void SetRootPreorderPartialsToStateFrequencies(const BeagleAccessories &ba) const;

  static inline void AddLowerPartialOperation(BeagleOperationVector &operations,
                                              const BeagleAccessories &ba, int node_id,
                                              int child0_id, int child1_id);
  static inline void AddUpperPartialOperation(BeagleOperationVector &operations,
                                              const BeagleAccessories &ba, int node_id,
                                              int sister_id, int parent_id);
  static inline std::pair<double, double> ComputeGradientEntry(
      BeagleAccessories &ba, const SizeVectorVector &indices_above, int node_id,
      int sister_id);
};

template <typename TOut, typename TTree, typename TTreeCollection>
std::vector<TOut> FatBeagleParallelize(
    std::function<TOut(FatBeagle *, const TTree &)> f,
    const std::vector<std::unique_ptr<FatBeagle>> &fat_beagles,
    const TTreeCollection &tree_collection, EigenMatrixXdRef param_matrix,
    const bool rescaling) {
  if (fat_beagles.empty()) {
    Failwith("Please add some FatBeagles that can be used for computation.");
  }
  std::vector<TOut> results(tree_collection.TreeCount());
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
      [&results, &tree_collection, &param_matrix, &rescaling, &f](FatBeagle *fat_beagle,
                                                                  size_t tree_number) {
        fat_beagle->SetParameters(param_matrix.row(tree_number));
        fat_beagle->SetRescaling(rescaling);
        results[tree_number] = f(fat_beagle, tree_collection.GetTree(tree_number));
      });
  return results;
}

// Tests live in unrooted_sbn_instance.hpp.
#endif  // SRC_FAT_BEAGLE_HPP_
