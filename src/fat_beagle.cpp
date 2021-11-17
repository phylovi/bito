// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "fat_beagle.hpp"

#include <numeric>
#include <utility>
#include <vector>

#include "rooted_gradient_transforms.hpp"

FatBeagle::FatBeagle(const PhyloModelSpecification &specification,
                     const SitePattern &site_pattern,
                     const FatBeagle::PackedBeagleFlags beagle_preference_flags,
                     bool use_tip_states)
    : phylo_model_(PhyloModel::OfSpecification(specification)),
      rescaling_(false),  // Note: rescaling_ set via the SetRescaling method.
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
double FatBeagle::LogLikelihoodInternals(
    const Node::NodePtr topology, const std::vector<double> &branch_lengths) const {
  BeagleAccessories ba(beagle_instance_, rescaling_, topology);
  BeagleOperationVector operations;
  beagleResetScaleFactors(beagle_instance_, 0);
  topology->BinaryIdPostorder(
      [&operations, &ba](int node_id, int child0_id, int child1_id) {
        AddLowerPartialOperation(operations, ba, node_id, child0_id, child1_id);
      });
  UpdateBeagleTransitionMatrices(ba, branch_lengths, nullptr);
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

double FatBeagle::LogLikelihood(const UnrootedTree &tree) const {
  auto detrifurcated_tree = tree.Detrifurcate();
  return LogLikelihoodInternals(detrifurcated_tree.Topology(),
                                detrifurcated_tree.BranchLengths());
}

double FatBeagle::UnrootedLogLikelihood(const RootedTree &tree) const {
  return LogLikelihoodInternals(tree.Topology(), tree.BranchLengths());
}

double LogDeterminantJacobian(const RootedTree &tree) {
  double log_det_jacobian = 0.0;
  size_t leaf_count = tree.LeafCount();
  tree.Topology()->TripleIdPreorderBifurcating(
      [&log_det_jacobian, &tree, leaf_count](size_t node_id, size_t sister_id,
                                             size_t parent_id) {
        if (node_id >= leaf_count) {
          log_det_jacobian +=
              std::log(tree.node_heights_[parent_id] - tree.node_bounds_[node_id]);
        }
      });
  return log_det_jacobian;
}

double FatBeagle::LogLikelihood(const RootedTree &tree) const {
  std::vector<double> branch_lengths = tree.BranchLengths();
  const std::vector<double> &rates = tree.GetRates();
  for (size_t i = 0; i < tree.BranchLengths().size() - 1; i++) {
    branch_lengths[i] *= rates[i];
  }
  return LogLikelihoodInternals(tree.Topology(), branch_lengths) +
         LogDeterminantJacobian(tree);
}

// Build differential matrix and scale it.
EigenMatrixXd BuildDifferentialMatrices(const SubstitutionModel &substitution_model,
                                        const EigenVectorXd &scalers) {
  size_t category_count = scalers.size();
  EigenMatrixXd Q = substitution_model.GetQMatrix();
  Eigen::Map<Eigen::RowVectorXd> mapQ(Q.data(), Q.size());
  EigenMatrixXd dQ = mapQ.replicate(category_count, 1);
  for (size_t k = 0; k < category_count; k++) {
    dQ.row(k) *= scalers[k];
  }
  return dQ;
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradientInternals(
    const Node::NodePtr topology, const std::vector<double> &branch_lengths,
    const EigenMatrixXd &dQ) const {
  beagleResetScaleFactors(beagle_instance_, 0);
  BeagleAccessories ba(beagle_instance_, rescaling_, topology);
  UpdateBeagleTransitionMatrices(ba, branch_lengths, nullptr);
  SetRootPreorderPartialsToStateFrequencies(ba);

  // Set differential matrices.
  int derivative_matrix_idx = ba.node_count_ - 1;
  beagleSetDifferentialMatrix(beagle_instance_, derivative_matrix_idx, dQ.data());
  const auto derivative_matrix_indices =
      std::vector<int>(ba.node_count_ - 1, derivative_matrix_idx);

  // Calculate post-order partials
  BeagleOperationVector operations;
  topology->BinaryIdPostorder(
      [&operations, &ba](int node_id, int child0_id, int child1_id) {
        AddLowerPartialOperation(operations, ba, node_id, child0_id, child1_id);
      });
  beagleUpdatePartials(beagle_instance_, operations.data(),
                       static_cast<int>(operations.size()),
                       ba.cumulative_scale_index_[0]);  // cumulative scale index

  // Calculate pre-order partials.
  operations.clear();
  topology->TripleIdPreorderBifurcating(
      [&operations, &ba](int node_id, int sister_id, int parent_id) {
        AddUpperPartialOperation(operations, ba, node_id, sister_id, parent_id);
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
      nullptr,                           // derivative-per-site output array
      gradient.data(),                   // sum of derivatives across sites output array
      nullptr);                          // sum of squared derivatives output array

  // Also calculate the likelihood.
  double log_like = 0.;
  beagleCalculateRootLogLikelihoods(
      beagle_instance_, &ba.root_id_, ba.category_weight_index_.data(),
      ba.state_frequency_index_.data(), ba.cumulative_scale_index_.data(),
      ba.mysterious_count_, &log_like);
  return {log_like, gradient};
}

FatBeagle *NullPtrAssert(FatBeagle *fat_beagle) {
  Assert(fat_beagle != nullptr, "NULL FatBeagle pointer!");
  return fat_beagle;
}

double FatBeagle::StaticUnrootedLogLikelihood(FatBeagle *fat_beagle,
                                              const UnrootedTree &in_tree) {
  return NullPtrAssert(fat_beagle)->LogLikelihood(in_tree);
}

double FatBeagle::StaticUnrootedLogLikelihoodOfRooted(FatBeagle *fat_beagle,
                                                      const RootedTree &in_tree) {
  return NullPtrAssert(fat_beagle)->UnrootedLogLikelihood(in_tree);
}

double FatBeagle::StaticRootedLogLikelihood(FatBeagle *fat_beagle,
                                            const RootedTree &in_tree) {
  return NullPtrAssert(fat_beagle)->LogLikelihood(in_tree);
}

PhyloGradient FatBeagle::StaticUnrootedGradient(FatBeagle *fat_beagle,
                                                const UnrootedTree &in_tree) {
  return NullPtrAssert(fat_beagle)->Gradient(in_tree);
}

PhyloGradient FatBeagle::StaticRootedGradient(FatBeagle *fat_beagle,
                                              const RootedTree &in_tree) {
  return NullPtrAssert(fat_beagle)->Gradient(in_tree);
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
  for (size_t i = 0; i < site_pattern.GetPatterns().size(); i++) {
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

void FatBeagle::UpdateSubstitutionModelInBeagle() const {
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

// If we pass nullptr as gradient_indices_ptr then we will not prepare for
// gradient calculation.
void FatBeagle::UpdateBeagleTransitionMatrices(
    const BeagleAccessories &ba, const std::vector<double> &branch_lengths,
    const int *const gradient_indices_ptr) const {
  beagleUpdateTransitionMatrices(beagle_instance_,         // instance
                                 0,                        // eigenIndex
                                 ba.node_indices_.data(),  // probabilityIndices
                                 gradient_indices_ptr,     // firstDerivativeIndices
                                 nullptr,                  // secondDerivativeIndices
                                 branch_lengths.data(),    // edgeLengths
                                 ba.node_count_ - 1);      // count
}

void FatBeagle::SetRootPreorderPartialsToStateFrequencies(
    const BeagleAccessories &ba) const {
  const EigenVectorXd &frequencies =
      phylo_model_->GetSubstitutionModel()->GetFrequencies();
  size_t category_count = phylo_model_->GetSiteModel()->GetCategoryCount();
  EigenVectorXd state_frequencies =
      frequencies.replicate(pattern_count_ * category_count, 1);
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
// Calculation of the substitution rate gradient.
// \partial{L}/\partial{r_i} = \partial{L}/\partial{b_i} \partial{b_i}/\partial{r_i}
// For strict clock:
// \partial{L}/\partial{r} = \sum_i \partial{L}/\partial{r_i}
std::vector<double> ClockGradient(const RootedTree &tree,
                                  const std::vector<double> &branch_gradient) {
  auto root_id = tree.Topology()->Id();
  std::vector<double> rate_gradient(root_id, 0);
  for (size_t i = 0; i < root_id; i++) {
    rate_gradient[i] = branch_gradient[i] * tree.branch_lengths_[i];
  }

  // Strict clock.
  if (tree.rate_count_ == 1) {
    return {std::accumulate(rate_gradient.cbegin(), rate_gradient.cend(), 0.0)};
  }
  // One rate per branch.
  else if (tree.rate_count_ == tree.rates_.size()) {
    return rate_gradient;
  } else {
    Failwith(
        "The number of rates should be equal to 1 (i.e. strict clock) or equal to "
        "the number of branches.");
  }
}

std::vector<double> DiscreteSiteModelGradient(
    const std::vector<double> &branch_lengths,
    const std::vector<double> &unscaled_category_gradient) {
  size_t edge_count = branch_lengths.size() - 1;
  double rate_gradient = 0;
  for (size_t node_id = 0; node_id < edge_count; node_id++) {
    rate_gradient += unscaled_category_gradient[node_id] * branch_lengths[node_id];
  }
  return {rate_gradient};
}

template <typename TTree>
std::vector<double> FatBeagle::SubstitutionModelGradientFiniteDifference(
    std::function<double(FatBeagle *, const TTree &)> f, FatBeagle *fat_beagle,
    const TTree &tree, SubstitutionModel *subst_model, const std::string &parameter_key,
    EigenVectorXd param_vector, double delta) const {
  return SubstitutionModelGradientFiniteDifference(f, fat_beagle, tree, subst_model,
                                                   parameter_key, param_vector, delta,
                                                   IdentityTransform());
}

template <typename TTree>
std::vector<double> FatBeagle::SubstitutionModelGradientFiniteDifference(
    std::function<double(FatBeagle *, const TTree &)> f, FatBeagle *fat_beagle,
    const TTree &tree, SubstitutionModel *subst_model, const std::string &parameter_key,
    EigenVectorXd param_vector, double delta, const Transform &transform) const {
  auto [parameter_start, parameter_length] =
      subst_model->GetBlockSpecification().GetMap().at(parameter_key);

  EigenVectorXd parameters = param_vector.segment(parameter_start, parameter_length);
  EigenVectorXd parameters_reparameterized = transform.inverse(parameters);

  std::vector<double> gradient(parameters_reparameterized.size());
  for (Eigen::Index parameter_idx = 0; parameter_idx < parameters_reparameterized.size();
       parameter_idx++) {
    double original_parameter_value = parameters_reparameterized[parameter_idx];
    parameters_reparameterized[parameter_idx] += delta;

    param_vector.segment(parameter_start, parameter_length) =
        transform(parameters_reparameterized);
    subst_model->SetParameters(param_vector);
    UpdateSubstitutionModelInBeagle();
    double log_prob_plus = f(fat_beagle, tree);

    parameters_reparameterized[parameter_idx] = original_parameter_value - delta;
    param_vector.segment(parameter_start, parameter_length) =
        transform(parameters_reparameterized);
    subst_model->SetParameters(param_vector);
    UpdateSubstitutionModelInBeagle();
    double log_prob_minus = f(fat_beagle, tree);

    gradient[parameter_idx] = (log_prob_plus - log_prob_minus) / (2. * delta);

    parameters_reparameterized[parameter_idx] = original_parameter_value;
    subst_model->SetParameters(param_vector);
  }
  UpdateSubstitutionModelInBeagle();
  return gradient;
}

template <typename TTree>
std::vector<double> FatBeagle::SubstitutionModelGradient(
    std::function<double(FatBeagle *, const TTree &)> f, FatBeagle *fat_beagle,
    const TTree &tree) const {
  auto subst_model = phylo_model_->GetSubstitutionModel();
  EigenVectorXd param_vector(subst_model->GetBlockSpecification().ParameterCount());
  auto subst_map = subst_model->GetBlockSpecification().GetMap();
  param_vector.segment(subst_map.at(SubstitutionModel::frequencies_key_).first,
                       subst_map.at(SubstitutionModel::frequencies_key_).second) =
      phylo_model_->GetSubstitutionModel()->GetFrequencies();
  param_vector.segment(subst_map.at(SubstitutionModel::rates_key_).first,
                       subst_map.at(SubstitutionModel::rates_key_).second) =
      phylo_model_->GetSubstitutionModel()->GetRates();
  // #324: make delta part of a gradient request
  double delta = 1.e-6;
  std::vector<double> frequencies_grad = SubstitutionModelGradientFiniteDifference(
      f, fat_beagle, tree, subst_model, SubstitutionModel::frequencies_key_,
      param_vector, delta, StickBreakingTransform());

  std::vector<double> gradient;
  // Rates in the GTR model are constrained to sum to 1
  if (subst_model->GetRates().size() == 6) {
    gradient = SubstitutionModelGradientFiniteDifference(
        f, fat_beagle, tree, subst_model, SubstitutionModel::rates_key_, param_vector,
        delta, StickBreakingTransform());
  } else {
    gradient = SubstitutionModelGradientFiniteDifference(
        f, fat_beagle, tree, subst_model, SubstitutionModel::rates_key_, param_vector,
        delta);
  }
  gradient.insert(gradient.end(), frequencies_grad.begin(), frequencies_grad.end());
  return gradient;
}

PhyloGradient FatBeagle::Gradient(const UnrootedTree &in_tree) const {
  auto tree = in_tree.Detrifurcate();
  tree.SlideRootPosition();
  EigenMatrixXd dQ =
      BuildDifferentialMatrices(*phylo_model_->GetSubstitutionModel(),
                                phylo_model_->GetSiteModel()->GetCategoryRates());
  auto [log_likelihood, branch_length_gradient] =
      BranchGradientInternals(tree.Topology(), tree.BranchLengths(), dQ);

  GradientMap gradient;

  // Calculate substitution model parameter gradient, if needed.
  if (phylo_model_->GetSubstitutionModel()->GetRates().size() > 0) {
    FatBeagle *mutable_this = const_cast<FatBeagle *>(this);
    gradient["substitution_model"] = SubstitutionModelGradient<UnrootedTree>(
        FatBeagle::StaticUnrootedLogLikelihood, mutable_this, in_tree);
  }

  auto site_model = phylo_model_->GetSiteModel();
  size_t category_count = site_model->GetCategoryCount();

  if (category_count > 1) {
    EigenMatrixXd dQ =
        BuildDifferentialMatrices(*phylo_model_->GetSubstitutionModel(),
                                  phylo_model_->GetSiteModel()->GetRateGradient());
    auto [log_likelihood, unscaled_category_gradient] =
        BranchGradientInternals(tree.Topology(), tree.BranchLengths(), dQ);
    gradient["site_model"] =
        DiscreteSiteModelGradient(tree.BranchLengths(), unscaled_category_gradient);
  }

  // We want the fixed node to have a zero gradient.
  branch_length_gradient[tree.Topology()->Children()[1]->Id()] = 0.;
  gradient["branch_lengths"] = branch_length_gradient;

  return {log_likelihood, gradient};
}

PhyloGradient FatBeagle::Gradient(const RootedTree &tree) const {
  // Scale time with clock rate.
  std::vector<double> branch_lengths = tree.BranchLengths();
  const std::vector<double> &rates = tree.GetRates();
  for (size_t i = 0; i < tree.BranchLengths().size() - 1; i++) {
    branch_lengths[i] *= rates[i];
  }

  // Calculate branch length gradient and log likelihood.
  EigenMatrixXd dQ =
      BuildDifferentialMatrices(*phylo_model_->GetSubstitutionModel(),
                                phylo_model_->GetSiteModel()->GetCategoryRates());
  auto [log_likelihood, branch_gradient] =
      BranchGradientInternals(tree.Topology(), branch_lengths, dQ);

  GradientMap gradient;

  gradient["branch_lengths"] = branch_gradient;
  // Calculate substitution model parameter gradient, if needed.
  if (phylo_model_->GetSubstitutionModel()->GetRates().size() > 0) {
    FatBeagle *mutable_this = const_cast<FatBeagle *>(this);
    gradient["substitution_model"] = SubstitutionModelGradient<RootedTree>(
        FatBeagle::StaticRootedLogLikelihood, mutable_this, tree);
  }

  // Calculate site model parameter gradient, if needed.
  auto site_model = phylo_model_->GetSiteModel();
  size_t category_count = site_model->GetCategoryCount();

  if (category_count > 1) {
    EigenMatrixXd dQ =
        BuildDifferentialMatrices(*phylo_model_->GetSubstitutionModel(),
                                  phylo_model_->GetSiteModel()->GetRateGradient());
    auto [log_likelihood, unscaled_category_gradient] =
        BranchGradientInternals(tree.Topology(), branch_lengths, dQ);
    gradient["site_model"] =
        DiscreteSiteModelGradient(branch_lengths, unscaled_category_gradient);
  }
  gradient["ratios_root_height"] = RatioGradientOfBranchGradient(tree, branch_gradient);
  gradient["clock_model"] = ClockGradient(tree, branch_gradient);

  return {log_likelihood, gradient};
}
