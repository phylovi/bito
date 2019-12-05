// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "fat_beagle.hpp"
#include <numeric>
#include <utility>
#include <vector>

FatBeagle::FatBeagle(const PhyloModelSpecification &specification,
                     const SitePattern &site_pattern)
    : phylo_model_(PhyloModel::OfSpecification(specification)),
      rescaling_(false),
      pattern_count_(static_cast<int>(site_pattern.PatternCount())) {
  beagle_instance_ = CreateInstance(site_pattern);
  SetTipStates(site_pattern);
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

double FatBeagle::LogLikelihood(const Tree &in_tree) const {
  beagleResetScaleFactors(beagle_instance_, 0);
  auto tree = PrepareTreeForLikelihood(in_tree);
  BeagleAccessories ba(beagle_instance_, rescaling_, tree);
  BeagleOperationVector operations;
  tree.Topology()->BinaryIdPostOrder(
      [&operations, &ba](int node_id, int child0_id, int child1_id) {
        AddLowerPartialOperation(operations, ba, node_id, child0_id, child1_id);
      });
  UpdateBeagleTransitionMatrices(ba, tree, nullptr);
  beagleUpdatePartials(beagle_instance_,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       ba.cumulative_scale_index_[0]);
  double log_like = 0;
  std::vector<int> root_id = {static_cast<int>(tree.Id())};
  beagleCalculateRootLogLikelihoods(
      beagle_instance_, root_id.data(), ba.category_weight_index_.data(),
      ba.state_frequency_index_.data(), ba.cumulative_scale_index_.data(),
      ba.mysterious_count_, &log_like);
  return log_like;
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradient(
    const Tree &in_tree) const {
  beagleResetScaleFactors(beagle_instance_, 0);
  auto tree = PrepareTreeForLikelihood(in_tree);
  tree.SlideRootPosition();

  BeagleAccessories ba(beagle_instance_, rescaling_, tree);
  BeagleOperationVector operations;
  const auto gradient_indices =
      BeagleAccessories::IotaVector(ba.node_count_ - 1, ba.node_count_);
  tree.Topology()->BinaryIdPostOrder(
      [&operations, &ba](int node_id, int child0_id, int child1_id) {
        AddLowerPartialOperation(operations, ba, node_id, child0_id, child1_id);
      });
  tree.Topology()->TripleIdPreOrderBifurcating(
      [&operations, &ba](int node_id, int sister_id, int parent_id) {
        AddUpperPartialOperation(operations, ba, node_id, sister_id, parent_id);
      });
  UpdateBeagleTransitionMatrices(ba, tree, gradient_indices.data());
  beagleUpdatePartials(beagle_instance_,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       BEAGLE_OP_NONE);  // cumulative scale index

  double log_like = 0;
  std::vector<double> gradient(ba.node_count_, 0.);
  SizeVectorVector indices_above = tree.Topology()->IdsAbove();
  for (auto &indices : indices_above) {
    // Reverse vector so we have indices from node to root.
    std::reverse(indices.begin(), indices.end());
    // Remove the root scalers.
    indices.pop_back();
  }
  // Actually compute gradient.
  tree.Topology()->TripleIdPreOrderBifurcating(
      [&ba, &gradient, &indices_above, &log_like](int node_id, int sister_id, int) {
        if (node_id != ba.fixed_node_id_) {
          auto [local_log_like, dlogLp] =
              ComputeGradientEntry(ba, indices_above, node_id, sister_id);
          log_like = local_log_like;
          gradient[static_cast<size_t>(node_id)] = dlogLp;
        }
      });
  return {log_like, gradient};
}

double FatBeagle::StaticLogLikelihood(FatBeagle *fat_beagle, const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "NULL FatBeagle pointer!");
  return fat_beagle->LogLikelihood(in_tree);
}

std::pair<double, std::vector<double>> FatBeagle::StaticBranchGradient(
    FatBeagle *fat_beagle, const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "NULL FatBeagle pointer!");
  return fat_beagle->BranchGradient(in_tree);
}

FatBeagle::BeagleInstance FatBeagle::CreateInstance(const SitePattern &site_pattern) {
  int taxon_count = static_cast<int>(site_pattern.SequenceCount());
  // Number of partial buffers to create (input):
  // taxon_count - 1 for lower partials (internal nodes only)
  // 2*taxon_count - 2 for upper partials (every node except the root)
  int partials_buffer_count = 3 * taxon_count - 3;
  // Number of compact state representation buffers to create -- for use with
  // setTipStates (input)
  int compact_buffer_count = taxon_count;
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
  int64_t preference_flags = 0;
  // Bit-flags indicating required implementation characteristics, see
  // BeagleFlags (input)
  int requirement_flags = BEAGLE_FLAG_SCALING_MANUAL;

  BeagleInstanceDetails return_info;
  auto beagle_instance = beagleCreateInstance(
      taxon_count, partials_buffer_count, compact_buffer_count, state_count,
      pattern_count, eigen_buffer_count, matrix_buffer_count, category_count,
      scale_buffer_count, allowed_resources, resource_count, preference_flags,
      requirement_flags, &return_info);
  if (return_info.flags & (BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_PROCESSOR_GPU)) {
    return beagle_instance;
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

Tree FatBeagle::PrepareTreeForLikelihood(const Tree &tree) {
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
  if (node_id != ba.root_child_id_ && node_id != ba.fixed_node_id_) {
    int upper_partial_index;
    int upper_matrix_index;
    if (parent_id == ba.root_child_id_) {
      upper_matrix_index = ba.root_child_id_;
      upper_partial_index = ba.fixed_node_id_;
    } else if (parent_id == ba.fixed_node_id_) {
      upper_matrix_index = ba.root_child_id_;
      upper_partial_index = ba.root_child_id_;
    } else {
      upper_partial_index = parent_id + ba.node_count_;
      upper_matrix_index = parent_id;
    }
    // Scalers are indexed differently for the upper conditional
    // likelihood. They start at the number of internal nodes + 1 because
    // of the lower conditional likelihoods. Also, in this case the leaves
    // have scalers.
    const int destinationScaleWrite =
        ba.rescaling_ ? node_id + 1 + ba.internal_count_ : BEAGLE_OP_NONE;
    operations.push_back({
        node_id + ba.node_count_,  // destinationPartials
        destinationScaleWrite, ba.destinationScaleRead_,
        upper_partial_index,  // child1Partials;
        upper_matrix_index,   // child1TransitionMatrix;
        sister_id,            // child2Partials;
        sister_id             // child2TransitionMatrix;
    });
  }
}

std::pair<double, double> FatBeagle::ComputeGradientEntry(
    BeagleAccessories &ba, const SizeVectorVector &indices_above, int node_id,
    int sister_id) {
  double log_like;
  double dlogLp;
  ba.upper_partials_index_[0] = node_id + ba.node_count_;
  ba.node_partial_indices_[0] = node_id;
  ba.node_mat_indices_[0] = node_id;
  ba.node_deriv_index_[0] = node_id + ba.node_count_;

  // If we're at the root then "up" is the sister.
  if (node_id == ba.root_child_id_) {
    ba.upper_partials_index_[0] = sister_id;
  }
  // Parent partial Buffers cannot be a taxon in
  // beagleCalculateEdgeLogLikelihoods.
  if (ba.node_partial_indices_[0] > ba.upper_partials_index_[0]) {
    std::swap(ba.node_partial_indices_, ba.upper_partials_index_);
  }

  if (ba.rescaling_) {
    beagleResetScaleFactors(ba.beagle_instance_, ba.cumulative_scale_index_[0]);
    auto scaler_indices =
        BeagleAccessories::IotaVector(static_cast<size_t>(ba.internal_count_ - 1), 1);
    // Replace lower scaler index with upper scaler index for nodes between
    // node_index and root.
    int child = node_id;
    for (size_t upper : indices_above[static_cast<size_t>(node_id)]) {
      int int_upper = static_cast<int>(upper);
      int scaler_indices_index = int_upper - ba.taxon_count_;
      Assert(scaler_indices_index >= 0, "int_upper must be >= taxon count.");
      scaler_indices[static_cast<size_t>(scaler_indices_index)] =
          child + ba.internal_count_ + 1;
      child = int_upper;
    }
    beagleAccumulateScaleFactors(ba.beagle_instance_, scaler_indices.data(),
                                 static_cast<int>(scaler_indices.size()),
                                 ba.cumulative_scale_index_[0]);
  }

  beagleCalculateEdgeLogLikelihoods(
      ba.beagle_instance_, ba.upper_partials_index_.data(),
      ba.node_partial_indices_.data(), ba.node_mat_indices_.data(),
      ba.node_deriv_index_.data(),
      nullptr,  // second derivative matrices
      ba.category_weight_index_.data(), ba.state_frequency_index_.data(),
      ba.cumulative_scale_index_.data(), ba.mysterious_count_, &log_like, &dlogLp,
      nullptr);  // destination for second derivative
  return {log_like, dlogLp};
}
