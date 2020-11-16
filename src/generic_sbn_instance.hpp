// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// This is a shared parent class for the rooted and unrooted SBN instances, templated on
// the type of trees and SBN support we need.
//
// The idea is to have as much in here as possible, which in practice means everything
// that doesn't require explicit reference to a specific tree type. Note that this
// excludes SBN training, which is different based on if we are thinking of trees as
// rooted or unrooted.

#ifndef SRC_GENERIC_SBN_INSTANCE_HPP_
#define SRC_GENERIC_SBN_INSTANCE_HPP_

#include "ProgressBar.hpp"
#include "alignment.hpp"
#include "csv.hpp"
#include "engine.hpp"
#include "mersenne_twister.hpp"
#include "numerical_utils.hpp"
#include "psp_indexer.hpp"
#include "rooted_sbn_support.hpp"
#include "sbn_probability.hpp"
#include "unrooted_sbn_support.hpp"

template <typename TTreeCollection, typename TSBNSupport,
          typename TIndexerRepresentation>
class GenericSBNInstance {
  static_assert(
      (std::is_same<TTreeCollection, RootedTreeCollection>::value &&
       std::is_same<TSBNSupport, RootedSBNSupport>::value &&
       std::is_same<TIndexerRepresentation, RootedIndexerRepresentation>::value) ||
      (std::is_same<TTreeCollection, UnrootedTreeCollection>::value &&
       std::is_same<TSBNSupport, UnrootedSBNSupport>::value &&
       std::is_same<TIndexerRepresentation, UnrootedIndexerRepresentation>::value));

 public:
  // A Range is a range of values of our vector using 0-indexed Python/Numpy slice
  // notation, such that if we have the range (1, 3), that refers to the items with
  // 0-index 1 and 2. Said another way, these are considered half-open intervals
  // [start, end).
  using Range = std::pair<size_t, size_t>;
  using RangeVector = std::vector<Range>;

  // The Primary Split Pair indexer.
  PSPIndexer psp_indexer_;
  // A vector that contains all of the SBN-related probabilities.
  EigenVectorXd sbn_parameters_;
  // The trees in our SBNInstance.
  TTreeCollection tree_collection_;

  // ** Initialization and status

  explicit GenericSBNInstance(std::string name)
      : name_(std::move(name)), rescaling_(false) {}

  size_t TaxonCount() const { return tree_collection_.TaxonCount(); }
  size_t TreeCount() const { return tree_collection_.TreeCount(); }
  const TagStringMap &TagTaxonMap() const { return tree_collection_.TagTaxonMap(); }
  const StringVector &TaxonNames() const { return sbn_support_.TaxonNames(); }
  const TSBNSupport &SBNSupport() const { return sbn_support_; }
  EigenConstVectorXdRef SBNParameters() const { return sbn_parameters_; }

  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (TreeCount()) {
      std::cout << TreeCount() << " unique tree topologies loaded on " << TaxonCount()
                << " leaves.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
    std::cout << alignment_.Data().size() << " sequences loaded.\n";
  }

  BitsetSizeDict RootsplitCounterOf(const Node::TopologyCounter &topologies) const {
    return TSBNSupport::RootsplitCounterOf(topologies);
  }
  PCSPCounter PCSPCounterOf(const Node::TopologyCounter &topologies) const {
    return TSBNSupport::PCSPCounterOf(topologies);
  }
  StringVector PrettyIndexer() const { return sbn_support_.PrettyIndexer(); }
  Node::TopologyCounter TopologyCounter() const {
    return tree_collection_.TopologyCounter();
  }

  // ** SBN-related items

  void SetSBNSupport(TSBNSupport &&sbn_support) {
    sbn_support_ = std::move(sbn_support);
    sbn_parameters_.resize(sbn_support_.GPCSPCount());
    sbn_parameters_.setOnes();
    psp_indexer_ = sbn_support_.BuildPSPIndexer();
  }

  // Use the loaded trees to set up the TopologyCounter, SBNSupport, etc.
  void ProcessLoadedTrees() {
    ClearTreeCollectionAssociatedState();
    topology_counter_ = TopologyCounter();
    SetSBNSupport(TSBNSupport(topology_counter_, tree_collection_.TaxonNames()));
  };

  // Any GPCSP that is not assigned a value by pretty_sbn_parameters will be assigned a
  // value of DOUBLE_MINIMUM (i.e. "log of 0"). We do emit a warning with this code is
  // used.
  //
  // We assume pretty_sbn_parameters is delivered in linear (i.e not log) space. If we
  // get log parameters they will have negative values, which will raise a failure.
  //
  // TODO am I happy with using a StringDoubleMap here?
  void SetSBNParameters(const StringDoubleMap &pretty_sbn_parameters,
                        bool warn_missing = true) {
    StringVector pretty_indexer = PrettyIndexer();
    size_t missing_count = 0;
    for (size_t i = 0; i < pretty_indexer.size(); i++) {  // NOLINT
      const std::string &pretty_gpcsp = pretty_indexer[i];
      const auto search = pretty_sbn_parameters.find(pretty_gpcsp);
      if (search == pretty_sbn_parameters.end()) {
        sbn_parameters_[i] = DOUBLE_MINIMUM;
        missing_count++;
      } else if (search->second > 0.) {
        sbn_parameters_[i] = log(search->second);
      } else if (search->second == 0.) {
        sbn_parameters_[i] = DOUBLE_MINIMUM;
      } else {
        Failwith(
            "Negative probability encountered in SetSBNParameters. Note that we don't "
            "expect the probabilities to be expressed in log space.");
      }
    }
    if (warn_missing && missing_count > 0) {
      std::cout << "Warning: when setting SBN parameters, " << missing_count
                << " were in the support but not specified; these were set to "
                << exp(DOUBLE_MINIMUM) << "." << std::endl;
    }
  }

  // The support has to already be set up to accept these SBN parameters.
  void ReadSBNParametersFromCSV(const std::string &csv_path) {
    SetSBNParameters(CSV::StringDoubleMapOfCSV(csv_path));
  }

  void CheckTopologyCounter() {
    if (TopologyCounter().empty()) {
      Failwith("Please load some trees into your SBN instance.");
    }
  }

  void CheckSBNSupportNonEmpty() {
    if (sbn_support_.Empty()) {
      Failwith("Please call ProcessLoadedTrees to prepare your SBN support.");
    }
  }

  void ProbabilityNormalizeSBNParametersInLog(EigenVectorXdRef sbn_parameters) const {
    sbn_support_.ProbabilityNormalizeSBNParametersInLog(sbn_parameters);
  }

  void TrainSimpleAverage() {
    CheckTopologyCounter();
    CheckSBNSupportNonEmpty();
    auto indexer_representation_counter =
        sbn_support_.IndexerRepresentationCounterOf(topology_counter_);
    SBNProbability::SimpleAverage(sbn_parameters_, indexer_representation_counter,
                                  sbn_support_.RootsplitCount(),
                                  sbn_support_.ParentToRange());
  }

  EigenVectorXd NormalizedSBNParameters() const {
    EigenVectorXd sbn_parameters_result = sbn_parameters_;
    ProbabilityNormalizeSBNParametersInLog(sbn_parameters_result);
    NumericalUtils::Exponentiate(sbn_parameters_result);
    return sbn_parameters_result;
  }

  StringDoubleVector PrettyIndexedSBNParameters() {
    StringDoubleVector result;
    auto sbn_parameters = NormalizedSBNParameters();
    result.reserve(sbn_parameters.size());
    const auto pretty_indexer = PrettyIndexer();
    for (size_t i = 0; i < pretty_indexer.size(); i++) {
      result.push_back({pretty_indexer.at(i), sbn_parameters(i)});
    }
    return result;
  }

  void SBNParametersToCSV(const std::string &file_path) {
    CSV::StringDoubleVectorToCSV(PrettyIndexedSBNParameters(), file_path);
  }

  // Get indexer representations of the trees in tree_collection_.
  // See the documentation of IndexerRepresentationOf in sbn_maps.hpp for an
  // explanation of what these are. This version uses the length of
  // sbn_parameters_ as a sentinel value for all rootsplits/PCSPs that aren't
  // present in the indexer.
  std::vector<TIndexerRepresentation> MakeIndexerRepresentations() const {
    std::vector<TIndexerRepresentation> representations;
    representations.reserve(tree_collection_.trees_.size());
    for (const auto &tree : tree_collection_.trees_) {
      representations.push_back(sbn_support_.IndexerRepresentationOf(tree.Topology()));
    }
    return representations;
  }

  // Calculate SBN probabilities for all currently-loaded trees.
  EigenVectorXd CalculateSBNProbabilities() {
    EigenVectorXd sbn_parameters_copy = sbn_parameters_;
    SBNProbability::ProbabilityNormalizeParamsInLog(sbn_parameters_copy,
                                                    sbn_support_.RootsplitCount(),
                                                    sbn_support_.ParentToRange());
    return SBNProbability::ProbabilityOfCollection(sbn_parameters_copy,
                                                   MakeIndexerRepresentations());
  }

  // ** Phylogenetic likelihood

  // Get the phylogenetic model parameters as a big matrix.
  Eigen::Ref<EigenMatrixXd> GetPhyloModelParams() { return phylo_model_params_; }

  // The phylogenetic model parameters broken down into blocks according to
  // model structure. See test_libsbn.py for an example of what this does.
  BlockSpecification::ParameterBlockMap GetPhyloModelParamBlockMap() {
    return GetEngine()->GetPhyloModelBlockSpecification().ParameterBlockMapOf(
        phylo_model_params_);
  }

  // Set whether we use rescaling for phylogenetic likelihood computation.
  void SetRescaling(bool use_rescaling) { rescaling_ = use_rescaling; }

  void CheckSequencesAndTreesLoaded() const {
    if (alignment_.SequenceCount() == 0) {
      Failwith(
          "Load an alignment into your SBNInstance on which you wish to "
          "calculate phylogenetic likelihoods.");
    }
    if (TreeCount() == 0) {
      Failwith(
          "Load some trees into your SBNInstance on which you wish to "
          "calculate phylogenetic likelihoods.");
    }
  }
  // Prepare for phylogenetic likelihood calculation. If we get a nullopt
  // argument, it just uses the number of trees currently in the SBNInstance.
  void PrepareForPhyloLikelihood(
      const PhyloModelSpecification &model_specification, size_t thread_count,
      const std::vector<BeagleFlags> &beagle_flag_vector = {},
      bool use_tip_states = true,
      const std::optional<size_t> &tree_count_option = std::nullopt) {
    const EngineSpecification engine_specification{thread_count, beagle_flag_vector,
                                                   use_tip_states};
    MakeEngine(engine_specification, model_specification);
    ResizePhyloModelParams(tree_count_option);
  }

  // Make the number of phylogentic model parameters fit the number of trees and
  // the speficied model. If we get a nullopt argument, it just uses the number
  // of trees currently in the SBNInstance.
  void ResizePhyloModelParams(std::optional<size_t> tree_count_option) {
    size_t tree_count = tree_count_option ? *tree_count_option : TreeCount();
    if (tree_count == 0) {
      Failwith(
          "Please add trees to your instance by sampling or loading before "
          "preparing for phylogenetic likelihood calculation.");
    }
    phylo_model_params_.resize(
        tree_count, GetEngine()->GetPhyloModelBlockSpecification().ParameterCount());
  }

  // ** I/O

  void ReadFastaFile(const std::string &fname) {
    alignment_ = Alignment::ReadFasta(fname);
  }

  // Allow users to pass in alignment directly.
  void SetAlignment(const Alignment &alignment) { alignment_ = alignment; }
  void SetAlignment(Alignment &&alignment) { alignment_ = alignment; }

 protected:
  // The name of our libsbn instance.
  std::string name_;
  // Our phylogenetic likelihood computation engine.
  std::unique_ptr<Engine> engine_;
  // Whether we use likelihood vector rescaling.
  bool rescaling_;
  // The multiple sequence alignment.
  Alignment alignment_;
  // The phylogenetic model parameterization. This has as many rows as there are
  // trees, and holds the parameters before likelihood computation, where they
  // will be processed across threads.
  EigenMatrixXd phylo_model_params_;
  // A counter for the currently loaded set of topologies.
  Node::TopologyCounter topology_counter_;
  TSBNSupport sbn_support_;

  MersenneTwister mersenne_twister_;
  inline void SetSeed(uint64_t seed) { mersenne_twister_.SetSeed(seed); }

  // Make a likelihood engine with the given specification.
  void MakeEngine(const EngineSpecification &engine_specification,
                  const PhyloModelSpecification &model_specification) {
    CheckSequencesAndTreesLoaded();
    SitePattern site_pattern(alignment_, TagTaxonMap());
    engine_ = std::make_unique<Engine>(engine_specification, model_specification,
                                       site_pattern);
  }

  // Return a raw pointer to the engine if it's available.
  Engine *GetEngine() const {
    if (engine_ != nullptr) {
      return engine_.get();
    }
    // else
    Failwith(
        "Engine not available. Call PrepareForPhyloLikelihood to make an "
        "engine for phylogenetic likelihood computation computation.");
  }

  // Sample an integer index in [range.first, range.second) according to
  // sbn_parameters_.
  size_t SampleIndex(Range range) const {
    const auto &[start, end] = range;
    Assert(start < end && end <= sbn_parameters_.size(),
           "SampleIndex given an invalid range.");
    // We do not want to overwrite sbn_parameters so we make a copy.
    EigenVectorXd sbn_parameters_subrange = sbn_parameters_.segment(start, end - start);
    NumericalUtils::ProbabilityNormalizeInLog(sbn_parameters_subrange);
    NumericalUtils::Exponentiate(sbn_parameters_subrange);
    std::discrete_distribution<> distribution(sbn_parameters_subrange.begin(),
                                              sbn_parameters_subrange.end());
    // We have to add on range.first because we have taken a slice of the full
    // array, and the sampler treats the beginning of this slice as zero.
    auto result =
        start + static_cast<size_t>(distribution(mersenne_twister_.GetGenerator()));
    Assert(result < end, "SampleIndex sampled a value out of range.");
    return result;
  }

  Node::NodePtr SampleTopology(bool rooted) const {
    // Start by sampling a rootsplit.
    size_t rootsplit_index =
        SampleIndex(std::pair<size_t, size_t>(0, sbn_support_.RootsplitCount()));
    const Bitset &rootsplit = sbn_support_.RootsplitsAt(rootsplit_index);
    // The addition below turns the rootsplit into a subsplit.
    auto topology = rooted ? SampleTopology(rootsplit + ~rootsplit)
                           : SampleTopology(rootsplit + ~rootsplit)->Deroot();
    topology->Polish();
    return topology;
  }

  // The input to this function is a parent subsplit (of length 2n).
  Node::NodePtr SampleTopology(const Bitset &parent_subsplit) const {
    auto process_subsplit = [this](const Bitset &parent) {
      auto singleton_option = parent.SplitChunk(1).SingletonOption();
      if (singleton_option) {
        return Node::Leaf(*singleton_option);
      }  // else
      auto child_index = SampleIndex(sbn_support_.ParentToRangeAt(parent));
      return SampleTopology(sbn_support_.IndexToChildAt(child_index));
    };
    return Node::Join(process_subsplit(parent_subsplit),
                      process_subsplit(parent_subsplit.RotateSubsplit()));
  }

  // Clear all of the state that depends on the current tree collection.
  void ClearTreeCollectionAssociatedState() {
    sbn_parameters_.resize(0);
    topology_counter_.clear();
    sbn_support_ = TSBNSupport();
  }

  void PushBackRangeForParentIfAvailable(const Bitset &parent,
                                         RangeVector &range_vector) {
    if (sbn_support_.ParentInSupport(parent)) {
      range_vector.push_back(sbn_support_.ParentToRangeAt(parent));
    }
  }

  RangeVector GetSubsplitRanges(const SizeVector &rooted_representation) {
    RangeVector subsplit_ranges;
    // PROFILE: should we be reserving here?
    subsplit_ranges.emplace_back(0, sbn_support_.RootsplitCount());
    Bitset root = sbn_support_.RootsplitsAt(rooted_representation[0]);
    PushBackRangeForParentIfAvailable(root + ~root, subsplit_ranges);
    PushBackRangeForParentIfAvailable(~root + root, subsplit_ranges);
    // Starting at 1 here because we took care of the rootsplit above (the 0th element).
    for (size_t i = 1; i < rooted_representation.size(); i++) {
      Bitset child = sbn_support_.IndexToChildAt(rooted_representation[i]);
      PushBackRangeForParentIfAvailable(child, subsplit_ranges);
      PushBackRangeForParentIfAvailable(child.RotateSubsplit(), subsplit_ranges);
    }
    return subsplit_ranges;
  }

  static EigenVectorXd CalculateMultiplicativeFactors(EigenVectorXdRef log_f) {
    double tree_count = log_f.size();
    double log_F = NumericalUtils::LogSum(log_f);
    double hat_L = log_F - log(tree_count);
    EigenVectorXd tilde_w = log_f.array() - log_F;
    tilde_w = tilde_w.array().exp();
    return hat_L - tilde_w.array();
  }

  static EigenVectorXd CalculateVIMCOMultiplicativeFactors(EigenVectorXdRef log_f) {
    // Use the geometric mean as \hat{f}(\tau^{-j}, \theta^{-j}), in eq:f_hat in
    // the implementation notes.
    size_t tree_count = log_f.size();
    double log_tree_count = log(tree_count);
    double sum_of_log_f = log_f.sum();
    // This has jth entry \hat{f}_{\bm{\phi},{\bm{\psi}}}(\tau^{-j},\bm{\theta}^{-j}),
    // i.e. the log of the geometric mean of each item other than j.
    EigenVectorXd log_geometric_mean =
        (sum_of_log_f - log_f.array()) / (tree_count - 1);
    EigenVectorXd per_sample_signal(tree_count);
    // This is a vector of entries that when summed become the parenthetical expression
    // in eq:perSampleLearning.
    EigenVectorXd log_f_perturbed = log_f;
    for (size_t j = 0; j < tree_count; j++) {
      log_f_perturbed(j) = log_geometric_mean(j);
      per_sample_signal(j) =
          log_f_perturbed.redux(NumericalUtils::LogAdd) - log_tree_count;
      // Reset the value.
      log_f_perturbed(j) = log_f(j);
    }
    EigenVectorXd multiplicative_factors = CalculateMultiplicativeFactors(log_f);
    multiplicative_factors -= per_sample_signal;
    return multiplicative_factors;
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED

#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_GENERIC_SBN_INSTANCE_HPP_
