// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <memory>
#include <queue>
#include <utility>
#include <vector>

#include "beagle_accessories.hpp"
#include "phylo_flags.hpp"
#include "phylo_model.hpp"
#include "phylo_gradient.hpp"
#include "rooted_tree_collection.hpp"
#include "site_pattern.hpp"
#include "stick_breaking_transform.hpp"
#include "task_processor.hpp"
#include "unrooted_tree_collection.hpp"

class FatBeagle {
 public:
  using PackedBeagleFlags = long;

  // This constructor makes the beagle_instance_
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

  double LogLikelihood(const UnrootedTree &tree,
                       std::optional<PhyloFlags> flags = std::nullopt) const;
  // This override performs a "classical" log likelihood calculation of a rooted tree
  // considered as an unrooted tree with no time-tree extras.
  double UnrootedLogLikelihood(const RootedTree &tree,
                               std::optional<PhyloFlags> flags = std::nullopt) const;
  double LogLikelihood(const RootedTree &tree,
                       std::optional<PhyloFlags> flags = std::nullopt) const;
  // Compute first derivative of the log likelihood with respect to each branch
  // length, as a vector of first derivatives indexed by node id.
  PhyloGradient Gradient(const UnrootedTree &tree,
                         std::optional<PhyloFlags> flags = std::nullopt) const;
  PhyloGradient Gradient(const RootedTree &tree,
                         std::optional<PhyloFlags> flags = std::nullopt) const;

  // ** Static Methods:
  // We can pass these static methods to FatBeagleParallelize.

  static double StaticUnrootedLogLikelihood(
      const FatBeagle *fat_beagle, const UnrootedTree &in_tree,
      std::optional<PhyloFlags> flags = std::nullopt);
  // This override performs a "classical" log likelihood calculation of a rooted tree
  // considered as an unrooted tree with no time-tree extras.
  static double StaticUnrootedLogLikelihoodOfRooted(
      const FatBeagle *fat_beagle, const RootedTree &in_tree,
      std::optional<PhyloFlags> flags = std::nullopt);
  static double StaticRootedLogLikelihood(
      const FatBeagle *fat_beagle, const RootedTree &in_tree,
      std::optional<PhyloFlags> flags = std::nullopt);
  static double StaticLogDetJacobianHeightTransform(
      const FatBeagle *fat_beagle, const RootedTree &in_tree,
      std::optional<PhyloFlags> flags = std::nullopt);

  static PhyloGradient StaticUnrootedGradient(
      const FatBeagle *fat_beagle, const UnrootedTree &in_tree,
      std::optional<PhyloFlags> flags = std::nullopt);
  static PhyloGradient StaticRootedGradient(
      const FatBeagle *fat_beagle, const RootedTree &in_tree,
      std::optional<PhyloFlags> flags = std::nullopt);
  static DoubleVector StaticGradientLogDeterminantJacobian(
      const FatBeagle *fat_beagle, const RootedTree &in_tree,
      std::optional<PhyloFlags> flags = std::nullopt);

  template <class TOut, class TTree>
  using StaticTreeFunction =
      std::function<TOut(const FatBeagle *, const TTree &, std::optional<PhyloFlags>)>;

  template <typename TTree>
  std::vector<double> SubstitutionModelGradientFiniteDifference(
      StaticTreeFunction<double, TTree> f, const FatBeagle *fat_beagle,
      const TTree &tree, SubstitutionModel *subst_model,
      const std::string &parameter_key, EigenVectorXd param_vector, double delta,
      std::optional<PhyloFlags> flags = std::nullopt) const;

  template <typename TTree>
  std::vector<double> SubstitutionModelGradientFiniteDifference(
      StaticTreeFunction<double, TTree> f, const FatBeagle *fat_beagle,
      const TTree &tree, SubstitutionModel *subst_model,
      const std::string &parameter_key, EigenVectorXd param_vector, double delta,
      const Transform &transform, std::optional<PhyloFlags> flags = std::nullopt) const;

  template <typename TTree>
  DoubleVectorPair SubstitutionModelGradient(
      StaticTreeFunction<double, TTree> f, const FatBeagle *fat_beagle,
      const TTree &tree, std::optional<PhyloFlags> flags = std::nullopt) const;

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
  void UpdateSubstitutionModelInBeagle() const;
  void UpdatePhyloModelInBeagle();

  double LogLikelihoodInternals(const Node::NodePtr topology,
                                const std::vector<double> &branch_lengths) const;
  std::pair<double, std::vector<double>> BranchGradientInternals(
      const Node::NodePtr topology, const std::vector<double> &branch_lengths,
      const EigenMatrixXd &dQ) const;

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
    FatBeagle::StaticTreeFunction<TOut, TTree> f,
    const std::vector<std::unique_ptr<FatBeagle>> &fat_beagles,
    const TTreeCollection &tree_collection, EigenMatrixXdRef param_matrix,
    const bool rescaling, std::optional<PhyloFlags> flags = std::nullopt) {
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
  Assert(static_cast<Eigen::Index>(tree_collection.TreeCount()) == param_matrix.rows(),
         "We param_matrix needs as many rows as we have trees.");

  TaskProcessor<FatBeagle *, size_t>(
      std::move(fat_beagle_queue), std::move(tree_number_queue),
      [&results, &tree_collection, &param_matrix, &rescaling, &f, &flags](
          FatBeagle *fat_beagle, size_t tree_number) {
        fat_beagle->SetParameters(param_matrix.row(tree_number));
        fat_beagle->SetRescaling(rescaling);
        results[tree_number] =
            f(fat_beagle, tree_collection.GetTree(tree_number), flags);
      });

  return results;
}

// Tests live in rooted_sbn_instance.hpp and unrooted_sbn_instance.hpp.
