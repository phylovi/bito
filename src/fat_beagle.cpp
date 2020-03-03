// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "fat_beagle.hpp"
#include <numeric>
#include <utility>
#include <vector>

FatBeagle::FatBeagle(const PhyloModelSpecification &specification,
                     const SitePattern &site_pattern,
                     const FatBeagle::PackedBeagleFlags beagle_preference_flags,
                     bool use_tip_states)
    : phylo_model_(PhyloModel::OfSpecification(specification)),
      rescaling_(false),
      pattern_count_(static_cast<int>(site_pattern.PatternCount())),
      use_tip_states_(use_tip_states) {
  std::tie(beagle_instance_, beagle_flags_) =
      CreateInstance(site_pattern, beagle_preference_flags);
  if (use_tip_states_) {
    SetTipStates(site_pattern);
  } else {
    SetTipPartials(site_pattern);
  }
  UpdatePhyloModelInBeagle();
};

FatBeagle::~FatBeagle() {
  auto finalize_result = beagleFinalizeInstance(beagle_instance_);
  if (finalize_result != 0) {
    std::cout << "beagleFinalizeInstance gave nonzero return value!";
    std::terminate();
  }
}

const BlockSpecification &FatBeagle::GetPhyloModelBlockSpecification() const {
  return phylo_model_->GetBlockSpecification();
}

void FatBeagle::SetParameters(const EigenVectorXdRef param_vector) {
  phylo_model_->SetParameters(param_vector);
  UpdatePhyloModelInBeagle();
}

// This is the "core" of the likelihood calculation, assuming that the tree is
// bifurcating.
double FatBeagle::LogLikelihoodInternals(const Tree &tree) const {
  BeagleAccessories ba(beagle_instance_, rescaling_, tree);
  BeagleOperationVector operations;
  beagleResetScaleFactors(beagle_instance_, 0);
  tree.Topology()->BinaryIdPostOrder(
      [&operations, &ba](int node_id, int child0_id, int child1_id) {
        AddLowerPartialOperation(operations, ba, node_id, child0_id, child1_id);
      });
  UpdateBeagleTransitionMatrices(ba, tree, nullptr);
  beagleUpdatePartials(beagle_instance_,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       ba.cumulative_scale_index_[0]);
  double log_like = 0.;
  beagleCalculateRootLogLikelihoods(
      beagle_instance_, &ba.root_id_, ba.category_weight_index_.data(),
      ba.state_frequency_index_.data(), ba.cumulative_scale_index_.data(),
      ba.mysterious_count_, &log_like);
  return log_like;
}

double FatBeagle::LogLikelihood(const Tree &tree) const {
  return LogLikelihoodInternals(DetrifurcateIfNeeded(tree));
}

double FatBeagle::LogLikelihood(const RootedTree &tree) const {
  return LogLikelihoodInternals(tree);
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradientInternals(
    const Tree &tree) const {
  beagleResetScaleFactors(beagle_instance_, 0);
  BeagleAccessories ba(beagle_instance_, rescaling_, tree);
  UpdateBeagleTransitionMatrices(ba, tree, nullptr);
  SetRootPreorderPartialsToStateFrequencies(ba);

  // Set differential matrix for each branch.
  const EigenMatrixXd &Q = phylo_model_->GetSubstitutionModel()->GetQMatrix();
  int derivative_matrix_idx = ba.node_count_ - 1;
  beagleSetDifferentialMatrix(beagle_instance_, derivative_matrix_idx, Q.data());
  const auto derivative_matrix_indices =
      std::vector<int>(ba.node_count_ - 1, derivative_matrix_idx);

  // Calculate post-order partials
  BeagleOperationVector operations;
  tree.Topology()->BinaryIdPostOrder(
      [&operations, &ba](int node_id, int child0_id, int child1_id) {
        AddLowerPartialOperation(operations, ba, node_id, child0_id, child1_id);
      });
  beagleUpdatePartials(beagle_instance_, operations.data(),
                       static_cast<int>(operations.size()),
                       ba.cumulative_scale_index_[0]);  // cumulative scale index

  // Calculate pre-order partials.
  operations.clear();
  tree.Topology()->TripleIdPreOrderBifurcating(
      [&operations, &ba](int node_id, int sister_id, int parent_id) {
        if (node_id != ba.root_id_) {
          AddUpperPartialOperation(operations, ba, node_id, sister_id, parent_id);
        }
      });
  beagleUpdatePrePartials(beagle_instance_, operations.data(),
                          static_cast<int>(operations.size()),
                          BEAGLE_OP_NONE);  // cumulative scale index

  // Actually compute the gradient.
  std::vector<double> gradient(ba.node_count_, 0.);
  const auto pre_buffer_indices =
      BeagleAccessories::IotaVector(ba.node_count_ - 1, ba.node_count_);
  beagleCalculateEdgeDerivatives(
      beagle_instance_,
      ba.node_indices_.data(),           // list of post order buffer indices
      pre_buffer_indices.data(),         // list of pre order buffer indices
      derivative_matrix_indices.data(),  // differential Q matrix indices
      ba.category_weight_index_.data(),  // category weights indices
      ba.node_count_ - 1,                // number of edges
      NULL,                              // derivative-per-site output array
      gradient.data(),                   // sum of derivatives across sites output array
      NULL);                             // sum of squared derivatives output array
  // We want the fixed node to have a zero gradient.
  gradient[ba.fixed_node_id_] = 0.;

  // Also calculate the likelihood.
  double log_like = 0.;
  beagleCalculateRootLogLikelihoods(
      beagle_instance_, &ba.root_id_, ba.category_weight_index_.data(),
      ba.state_frequency_index_.data(), ba.cumulative_scale_index_.data(),
      ba.mysterious_count_, &log_like);
  return {log_like, gradient};
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradient(
    const Tree &in_tree) const {
  auto tree = DetrifurcateIfNeeded(in_tree);
  tree.SlideRootPosition();
  return BranchGradientInternals(tree);
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradient(
    const RootedTree &in_tree) const {
  std::cout << "We are doing a ROOTED branch gradient calculation.\n";
  return BranchGradientInternals(in_tree);
}

FatBeagle *NullPtrAssert(FatBeagle *fat_beagle) {
  Assert(fat_beagle != nullptr, "NULL FatBeagle pointer!");
  return fat_beagle;
}

double FatBeagle::StaticLogLikelihood(FatBeagle *fat_beagle, const Tree &in_tree) {
  return NullPtrAssert(fat_beagle)->LogLikelihood(in_tree);
}

double FatBeagle::StaticRootedLogLikelihood(FatBeagle *fat_beagle,
                                            const RootedTree &in_tree) {
  return NullPtrAssert(fat_beagle)->LogLikelihood(in_tree);
}

std::pair<double, std::vector<double>> FatBeagle::StaticBranchGradient(
    FatBeagle *fat_beagle, const Tree &in_tree) {
  return NullPtrAssert(fat_beagle)->BranchGradient(in_tree);
}

std::pair<double, std::vector<double>> FatBeagle::StaticRootedBranchGradient(
    FatBeagle *fat_beagle, const RootedTree &in_tree) {
  return NullPtrAssert(fat_beagle)->BranchGradient(in_tree);
}

std::pair<FatBeagle::BeagleInstance, FatBeagle::PackedBeagleFlags>
FatBeagle::CreateInstance(const SitePattern &site_pattern,
                          FatBeagle::PackedBeagleFlags beagle_preference_flags) {
  int taxon_count = static_cast<int>(site_pattern.SequenceCount());
  // Number of partial buffers to create (input):
  // taxon_count - 1 for lower partials (internal nodes only)
  // 2*taxon_count - 1 for upper partials (every node)
  int partials_buffer_count = 3 * taxon_count - 2;
  if (!use_tip_states_) {
    partials_buffer_count += taxon_count;
  }
  // Number of compact state representation buffers to create -- for use with
  // setTipStates (input)
  int compact_buffer_count = (use_tip_states_ ? taxon_count : 0);
  // The number of states.
  int state_count =
      static_cast<int>(phylo_model_->GetSubstitutionModel()->GetStateCount());
  // Number of site patterns to be handled by the instance.
  int pattern_count = pattern_count_;
  // Number of eigen-decomposition buffers to allocate (input)
  int eigen_buffer_count = 1;
  // Number of transition matrix buffers (input) -- two per edge
  int matrix_buffer_count = 2 * (2 * taxon_count - 1);
  // Number of rate categories
  int category_count =
      static_cast<int>(phylo_model_->GetSiteModel()->GetCategoryCount());
  // Number of scaling buffers -- 1 buffer per partial buffer and 1 more
  // for accumulating scale factors in position 0.
  int scale_buffer_count = partials_buffer_count + 1;
  // List of potential resources on which this instance is allowed (input,
  // NULL implies no restriction
  int *allowed_resources = nullptr;
  // Length of resourceList list (input) -- not needed to use the default
  // hardware config
  int resource_count = 0;
  // Bit-flags indicating preferred implementation charactertistics, see
  // BeagleFlags (input)
  int requirement_flags = BEAGLE_FLAG_SCALING_MANUAL;

  BeagleInstanceDetails return_info;
  auto beagle_instance = beagleCreateInstance(
      taxon_count, partials_buffer_count, compact_buffer_count, state_count,
      pattern_count, eigen_buffer_count, matrix_buffer_count, category_count,
      scale_buffer_count, allowed_resources, resource_count, beagle_preference_flags,
      requirement_flags, &return_info);
  if (return_info.flags & (BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_PROCESSOR_GPU)) {
    return {beagle_instance, return_info.flags};
  }  // else
  Failwith("Couldn't get a CPU or a GPU from BEAGLE.");
}

void FatBeagle::SetTipStates(const SitePattern &site_pattern) {
  int taxon_number = 0;
  for (const auto &pattern : site_pattern.GetPatterns()) {
    beagleSetTipStates(beagle_instance_, taxon_number++, pattern.data());
  }
  beagleSetPatternWeights(beagle_instance_, site_pattern.GetWeights().data());
}

void FatBeagle::SetTipPartials(const SitePattern &site_pattern) {
  for (int i = 0; i < site_pattern.GetPatterns().size(); i++) {
    beagleSetTipPartials(beagle_instance_, i, site_pattern.GetPartials(i).data());
  }
  beagleSetPatternWeights(beagle_instance_, site_pattern.GetWeights().data());
}

void FatBeagle::UpdateSiteModelInBeagle() {
  const auto &site_model = phylo_model_->GetSiteModel();
  const auto &weights = site_model->GetCategoryProportions();
  const auto &rates = site_model->GetCategoryRates();
  beagleSetCategoryWeights(beagle_instance_, 0, weights.data());
  beagleSetCategoryRates(beagle_instance_, rates.data());
}

void FatBeagle::UpdateSubstitutionModelInBeagle() {
  const auto &substitution_model = phylo_model_->GetSubstitutionModel();
  const EigenMatrixXd &eigenvectors = substitution_model->GetEigenvectors();
  const EigenMatrixXd &inverse_eigenvectors =
      substitution_model->GetInverseEigenvectors();
  const EigenVectorXd &eigenvalues = substitution_model->GetEigenvalues();
  const EigenVectorXd &frequencies = substitution_model->GetFrequencies();

  beagleSetStateFrequencies(beagle_instance_, 0, frequencies.data());
  beagleSetEigenDecomposition(beagle_instance_,
                              0,  // eigenIndex
                              &eigenvectors.data()[0], &inverse_eigenvectors.data()[0],
                              &eigenvalues.data()[0]);
}

void FatBeagle::UpdatePhyloModelInBeagle() {
  // Issue #146: put in a clock model here.
  UpdateSiteModelInBeagle();
  UpdateSubstitutionModelInBeagle();
}

Tree FatBeagle::DetrifurcateIfNeeded(const Tree &tree) {
  if (tree.Children().size() == 3) {
    return tree.Detrifurcate();
  }  // else
  if (tree.Children().size() == 2) {
    return tree;
  }  // else
  Failwith(
      "Tree likelihood calculations should be done on a tree with a "
      "bifurcation or a trifurcation at the root.");
}

// If we pass nullptr as gradient_indices_ptr then we will not prepare for
// gradient calculation.
void FatBeagle::UpdateBeagleTransitionMatrices(
    const BeagleAccessories &ba, const Tree &tree,
    const int *const gradient_indices_ptr) const {
  beagleUpdateTransitionMatrices(beagle_instance_,         // instance
                                 0,                        // eigenIndex
                                 ba.node_indices_.data(),  // probabilityIndices
                                 gradient_indices_ptr,     // firstDerivativeIndices
                                 nullptr,                  // secondDerivativeIndices
                                 tree.BranchLengths().data(),  // edgeLengths
                                 ba.node_count_ - 1);          // count
}

void FatBeagle::SetRootPreorderPartialsToStateFrequencies(
    const BeagleAccessories &ba) const {
  const EigenVectorXd &frequencies =
      phylo_model_->GetSubstitutionModel()->GetFrequencies();
  EigenVectorXd state_frequencies = frequencies.replicate(pattern_count_, 1);
  beagleSetPartials(beagle_instance_, ba.root_id_ + ba.node_count_,
                    state_frequencies.data());
}

void FatBeagle::AddLowerPartialOperation(BeagleOperationVector &operations,
                                         const BeagleAccessories &ba, const int node_id,
                                         const int child0_id, const int child1_id) {
  const int destinationScaleWrite =
      ba.rescaling_ ? node_id - ba.taxon_count_ + 1 : BEAGLE_OP_NONE;
  // We can't emplace_back because BeagleOperation has no constructor.
  // The compiler should elide this though.
  operations.push_back({
      node_id,  // destinationPartials
      destinationScaleWrite, ba.destinationScaleRead_,
      child0_id,  // child1Partials;
      child0_id,  // child1TransitionMatrix;
      child1_id,  // child2Partials;
      child1_id   // child2TransitionMatrix;
  });
}

void FatBeagle::AddUpperPartialOperation(BeagleOperationVector &operations,
                                         const BeagleAccessories &ba, const int node_id,
                                         const int sister_id, const int parent_id) {
  // Scalers are indexed differently for the upper conditional
  // likelihood. They start at the number of internal nodes + 1 because
  // of the lower conditional likelihoods. Also, in this case the leaves
  // have scalers.
  const int destinationScaleWrite =
      ba.rescaling_ ? node_id + 1 + ba.internal_count_ : BEAGLE_OP_NONE;

  operations.push_back({
      node_id + ba.node_count_,  // dest pre-order partial of current node
      destinationScaleWrite, ba.destinationScaleRead_,
      parent_id + ba.node_count_,  // pre-order partial parent
      node_id,                     // matrices of current node
      sister_id,                   // post-order partial of sibling
      sister_id                    // matrices of sibling
  });
}
