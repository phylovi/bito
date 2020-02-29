// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"
#include <iostream>
#include <memory>
#include <unordered_set>

#include "eigen_sugar.hpp"
#include "numerical_utils.hpp"

void SBNInstance::PrintStatus() {
  std::cout << "Status for instance '" << name_ << "':\n";
  if (tree_collection_.TreeCount()) {
    std::cout << TreeCount() << " unique tree topologies loaded on "
              << tree_collection_.TaxonCount() << " leaves.\n";
  } else {
    std::cout << "No trees loaded.\n";
  }
  std::cout << alignment_.Data().size() << " sequences loaded.\n";
}

// ** Building SBN-related items

void SBNInstance::ProcessLoadedTrees() {
  size_t index = 0;
  ClearTreeCollectionAssociatedState();
  topology_counter_ = tree_collection_.TopologyCounter();
  // Start by adding the rootsplits.
  for (const auto &iter : SBNMaps::RootsplitCounterOf(topology_counter_)) {
    SafeInsert(indexer_, iter.first, index);
    rootsplits_.push_back(iter.first);
    index++;
  }
  // Now add the PCSSs.
  for (const auto &[parent, child_counter] :
       SBNMaps::PCSSCounterOf(topology_counter_)) {
    SafeInsert(parent_to_range_, parent, {index, index + child_counter.size()});
    for (const auto &child_iter : child_counter) {
      const auto &child = child_iter.first;
      SafeInsert(indexer_, parent + child, index);
      SafeInsert(index_to_child_, index, Bitset::ChildSubsplit(parent, child));
      index++;
    }
  }
  sbn_parameters_.resize(index);
  sbn_parameters_.setOnes();
  psp_indexer_ = PSPIndexer(rootsplits_, indexer_);
  taxon_names_ = tree_collection_.TaxonNames();
}

void SBNInstance::CheckSBNMapsAvailable() {
  if (indexer_.empty() || index_to_child_.empty() || parent_to_range_.empty() ||
      rootsplits_.empty() || taxon_names_.empty()) {
    Failwith("Please call ProcessLoadedTrees to prepare your SBN maps.");
  }
}

StringVector SBNInstance::PrettyIndexer() {
  StringVector pretty_representation(indexer_.size());
  for (const auto &[key, idx] : indexer_) {
    if (idx < rootsplits_.size()) {
      pretty_representation[idx] = key.ToString();
    } else {
      pretty_representation[idx] = key.PCSSToString();
    }
  }
  return pretty_representation;
}

void SBNInstance::PrettyPrintIndexer() {
  auto pretty_representation = PrettyIndexer();
  for (size_t i = 0; i < pretty_representation.size(); i++) {
    std::cout << i << "\t" << pretty_representation[i] << std::endl;
  }
}

void SBNInstance::TrainSimpleAverage() {
  auto indexer_representation_counter = SBNMaps::IndexerRepresentationCounterOf(
      indexer_, topology_counter_, sbn_parameters_.size());
  SBNProbability::SimpleAverage(sbn_parameters_, indexer_representation_counter,
                                rootsplits_.size(), parent_to_range_);
}

EigenVectorXd SBNInstance::TrainExpectationMaximization(double alpha, size_t max_iter,
                                                        double score_epsilon) {
  auto indexer_representation_counter = SBNMaps::IndexerRepresentationCounterOf(
      indexer_, topology_counter_, sbn_parameters_.size());
  return SBNProbability::ExpectationMaximization(
      sbn_parameters_, indexer_representation_counter, rootsplits_.size(),
      parent_to_range_, alpha, max_iter, score_epsilon);
}

EigenVectorXd SBNInstance::CalculateSBNProbabilities() {
  EigenVectorXd sbn_parameters_copy = sbn_parameters_;
  SBNProbability::ProbabilityNormalizeParamsInLog(sbn_parameters_copy,
                                                  rootsplits_.size(), parent_to_range_);
  return SBNProbability::ProbabilityOf(sbn_parameters_copy,
                                       MakeIndexerRepresentations());
}

size_t SBNInstance::SampleIndex(std::pair<size_t, size_t> range) const {
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
  auto result = start + static_cast<size_t>(distribution(random_generator_));
  Assert(result < end, "SampleIndex sampled a value out of range.");
  return result;
}

// This function samples a tree by first sampling the rootsplit, and then
// calling the recursive form of SampleTopology.
Node::NodePtr SBNInstance::SampleTopology(bool rooted) const {
  // Start by sampling a rootsplit.
  size_t rootsplit_index =
      SampleIndex(std::pair<size_t, size_t>(0, rootsplits_.size()));
  const Bitset &rootsplit = rootsplits_.at(rootsplit_index);
  // The addition below turns the rootsplit into a subsplit.
  auto topology = rooted ? SampleTopology(rootsplit + ~rootsplit)
                         : SampleTopology(rootsplit + ~rootsplit)->Deroot();
  topology->Polish();
  return topology;
}

// The input to this function is a parent subsplit (of length 2n).
Node::NodePtr SBNInstance::SampleTopology(const Bitset &parent_subsplit) const {
  auto process_subsplit = [this](const Bitset &parent) {
    auto singleton_option = parent.SplitChunk(1).SingletonOption();
    if (singleton_option) {
      return Node::Leaf(*singleton_option);
    }  // else
    auto child_index = SampleIndex(parent_to_range_.at(parent));
    return SampleTopology(index_to_child_.at(child_index));
  };
  return Node::Join(process_subsplit(parent_subsplit),
                    process_subsplit(parent_subsplit.RotateSubsplit()));
}

void SBNInstance::SampleTrees(size_t count) {
  CheckSBNMapsAvailable();
  auto leaf_count = rootsplits_[0].size();
  // 2n-2 because trees are unrooted.
  auto edge_count = 2 * static_cast<int>(leaf_count) - 2;
  tree_collection_.trees_.clear();
  for (size_t i = 0; i < count; i++) {
    std::vector<double> branch_lengths(static_cast<size_t>(edge_count));
    tree_collection_.trees_.emplace_back(
        Tree(SampleTopology(), std::move(branch_lengths)));
  }
}

std::vector<IndexerRepresentation> SBNInstance::MakeIndexerRepresentations() const {
  std::vector<IndexerRepresentation> representations;
  representations.reserve(tree_collection_.trees_.size());
  for (const auto &tree : tree_collection_.trees_) {
    representations.push_back(SBNMaps::IndexerRepresentationOf(
        indexer_, tree.Topology(), sbn_parameters_.size()));
  }
  return representations;
}

std::vector<SizeVectorVector> SBNInstance::MakePSPIndexerRepresentations() const {
  std::vector<SizeVectorVector> representations;
  representations.reserve(tree_collection_.trees_.size());
  for (const auto &tree : tree_collection_.trees_) {
    representations.push_back(psp_indexer_.RepresentationOf(tree.Topology()));
  }
  return representations;
}

StringVector SBNInstance::StringReversedIndexer() const {
  std::vector<std::string> reversed_indexer(indexer_.size());
  for (const auto &[key, idx] : indexer_) {
    if (idx < rootsplits_.size()) {
      reversed_indexer[idx] = key.ToString();
    } else {
      reversed_indexer[idx] = key.PCSSToString();
    }
  }
  return reversed_indexer;
}

StringSetVector SBNInstance::StringIndexerRepresentationOf(
    IndexerRepresentation indexer_representation) const {
  auto reversed_indexer = StringReversedIndexer();
  StringSetVector string_sets;
  for (const auto &rooted_representation : indexer_representation) {
    StringSet string_set;
    for (const auto index : rooted_representation) {
      SafeInsert(string_set, reversed_indexer[index]);
    }
    string_sets.push_back(std::move(string_set));
  }
  return string_sets;
}

DoubleVectorVector SBNInstance::SplitLengths() const {
  return psp_indexer_.SplitLengths(tree_collection_);
}

// ** I/O

std::tuple<StringSizeMap, StringSizePairMap> SBNInstance::GetIndexers() const {
  auto str_indexer = StringifyMap(indexer_);
  auto str_parent_to_range = StringifyMap(parent_to_range_);
  std::string rootsplit("rootsplit");
  SafeInsert(str_parent_to_range, rootsplit, {0, rootsplits_.size()});
  return {str_indexer, str_parent_to_range};
}

// This function is really just for testing-- it recomputes from scratch.
std::pair<StringSizeMap, StringPCSSMap> SBNInstance::SplitCounters() const {
  auto counter = tree_collection_.TopologyCounter();
  return {StringifyMap(SBNMaps::RootsplitCounterOf(counter).Map()),
          SBNMaps::StringPCSSMapOf(SBNMaps::PCSSCounterOf(counter))};
}

void SBNInstance::ReadNewickFile(std::string fname) {
  Driver driver;
  tree_collection_ = driver.ParseNewickFile(fname);
}

void SBNInstance::ReadNexusFile(std::string fname) {
  Driver driver;
  tree_collection_ = driver.ParseNexusFile(fname);
}

void SBNInstance::ReadFastaFile(std::string fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

// ** Phylogenetic likelihood

void SBNInstance::CheckSequencesAndTreesLoaded() const {
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

Eigen::Ref<EigenMatrixXd> SBNInstance::GetPhyloModelParams() {
  return phylo_model_params_;
}

BlockSpecification::ParameterBlockMap SBNInstance::GetPhyloModelParamBlockMap() {
  return GetEngine()->GetPhyloModelBlockSpecification().ParameterBlockMapOf(
      phylo_model_params_);
}

void SBNInstance::MakeEngine(const EngineSpecification &engine_specification,
                             const PhyloModelSpecification &model_specification) {
  CheckSequencesAndTreesLoaded();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
  engine_ =
      std::make_unique<Engine>(engine_specification, model_specification, site_pattern);
}

Engine *SBNInstance::GetEngine() const {
  if (engine_ != nullptr) {
    return engine_.get();
  }
  // else
  Failwith(
      "Engine not available. Call PrepareForPhyloLikelihood to make an "
      "engine for phylogenetic likelihood computation computation.");
}

void SBNInstance::ClearTreeCollectionAssociatedState() {
  sbn_parameters_.resize(0);
  rootsplits_.clear();
  indexer_.clear();
  index_to_child_.clear();
  parent_to_range_.clear();
  topology_counter_.clear();
}

void SBNInstance::PrepareForPhyloLikelihood(
    const PhyloModelSpecification &model_specification, size_t thread_count,
    const std::vector<BeagleFlags> &beagle_flag_vector, const bool use_tip_states,
    std::optional<size_t> tree_count_option) {
  const EngineSpecification engine_specification{thread_count, beagle_flag_vector,
                                                 use_tip_states};
  MakeEngine(engine_specification, model_specification);
  ResizePhyloModelParams(tree_count_option);
}

void SBNInstance::ResizePhyloModelParams(std::optional<size_t> tree_count_option) {
  size_t tree_count =
      tree_count_option ? *tree_count_option : tree_collection_.TreeCount();
  if (tree_count == 0) {
    Failwith(
        "Please add trees to your instance by sampling or loading before "
        "preparing for phylogenetic likelihood calculation.");
  }
  phylo_model_params_.resize(
      tree_count, GetEngine()->GetPhyloModelBlockSpecification().ParameterCount());
}

std::vector<double> SBNInstance::LogLikelihoods() {
  return GetEngine()->LogLikelihoods(tree_collection_, phylo_model_params_, rescaling_);
}

std::vector<std::pair<double, std::vector<double>>> SBNInstance::BranchGradients() {
  return GetEngine()->BranchGradients(tree_collection_, phylo_model_params_,
                                      rescaling_);
}

// Retrieves range of subsplits for each s|t that appears in the tree
// given by rooted_representation.
std::vector<std::pair<size_t, size_t>> SBNInstance::GetSubsplitRanges(
    const SizeVector &rooted_representation) {
  std::vector<std::pair<size_t, size_t>> subsplits;
  subsplits.push_back(std::make_pair(0, rootsplits_.size()));
  Bitset root = rootsplits_.at(rooted_representation[0]);
  // add child subsplit ranges
  if (parent_to_range_.count(root + ~root) > 0) {
    const auto &parent_range = parent_to_range_.at(root + ~root);
    subsplits.push_back(parent_range);
  }
  if (parent_to_range_.count(~root + root) > 0) {
    const auto &parent_range = parent_to_range_.at(~root + root);
    subsplits.push_back(parent_range);
  }
  for (size_t i = 1; i < rooted_representation.size(); i++) {
    size_t idx = rooted_representation[i];
    Bitset parent = index_to_child_.at(idx);
    if (parent_to_range_.count(parent) > 0) {
      const auto &parent_range = parent_to_range_.at(parent);
      subsplits.push_back(parent_range);
    }
    if (parent_to_range_.count(parent.RotateSubsplit()) > 0) {
      const auto &parent_range = parent_to_range_.at(parent.RotateSubsplit());
      subsplits.push_back(parent_range);
    }
  }
  return subsplits;
}

// This multiplicative factor is the quantity inside the parentheses in eq:nabla in the
// tex.
EigenVectorXd CalculateMultiplicativeFactors(const EigenVectorXdRef log_f) {
  size_t num_trees = log_f.size();
  double log_F = NumericalUtils::LogSum(log_f);
  double hat_L = log_F - log(num_trees);
  EigenVectorXd tilde_w = log_f.array() - log_F;
  tilde_w = tilde_w.array().exp();
  EigenVectorXd multiplicative_factors = hat_L - tilde_w.array();
  return multiplicative_factors;
}

// This gives the gradient of log q at a specific unrooted topology.
// See eq:gradLogQ in the tex.
EigenVectorXd SBNInstance::GradientOfLogQ(
                          const IndexerRepresentation &indexer_representation) {
  EigenVectorXd grad_log_q = EigenVectorXd::Zero(sbn_parameters_.size());
  double log_q = DOUBLE_NEG_INF;
  for (const auto &rooted_representation : indexer_representation) {
    if (!SBNProbability::IsInSBNSupport(rooted_representation,
                                        sbn_parameters_.size())) {
      continue;
    }
    double log_probability_rooted_tree =
        SBNProbability::SumOf(sbn_parameters_, rooted_representation, 0.0);
    // We need to look up the subsplits in the tree. Set representation allows
    // fast lookup.
    std::unordered_set<size_t> rooted_representation_as_set(
      rooted_representation.begin(), rooted_representation.end());
    // Get all subsplit ranges.
    auto subsplit_ranges = GetSubsplitRanges(rooted_representation);
    // Now, we update the gradients.
    for (const auto &[begin, end] : subsplit_ranges) {
      for (size_t idx = begin; idx < end; idx++) {
        size_t ind_subsplit_in_rooted_tree =
          (rooted_representation_as_set.count(idx) > 0) ? 1 : 0;
        // %EM What if, instead of this if, we actually had an indicator variable? I
        // feel like this might be easier to follow for people who are reading the tex.
        // %SJ Done.
        grad_log_q[idx] += exp(log_probability_rooted_tree) *
            (ind_subsplit_in_rooted_tree - exp(sbn_parameters_[idx]));
      }
    }
    log_q = NumericalUtils::LogAdd(log_q, log_probability_rooted_tree);
  }
  grad_log_q = grad_log_q.array() * exp(-log_q);
  return grad_log_q;
}

EigenVectorXd SBNInstance::TopologyGradients(const EigenVectorXdRef log_f) {
  size_t num_trees = tree_collection_.Trees().size();

  EigenVectorXd gradient_vector = EigenVectorXd::Zero(sbn_parameters_.size());

  EigenVectorXd multiplicative_factors = CalculateMultiplicativeFactors(log_f);

  // %EM I'm not sure what this comment means:
  // %SJ How about now?
  // Normalize the sbn parameters in log space.
  NormalizeSBNParametersInLog();
  for (size_t i = 0; i < num_trees; i++) {
    double multiplicative_factor = multiplicative_factors(i);
    const auto indexer_representation = SBNMaps::IndexerRepresentationOf(
      indexer_, tree_collection_.GetTree(i).Topology(), sbn_parameters_.size());
    EigenVectorXd log_grad_q = GradientOfLogQ(indexer_representation);
    log_grad_q = log_grad_q.array() * multiplicative_factor;
    gradient_vector = gradient_vector + log_grad_q;
  }
  return gradient_vector;
}

void SBNInstance::NormalizeSBNParametersInLog() {
  SBNProbability::ProbabilityNormalizeParamsInLog(sbn_parameters_,
    rootsplits_.size(), parent_to_range_);
}

// Here we initialize our static random number generator.
std::random_device SBNInstance::random_device_;
std::mt19937 SBNInstance::random_generator_(SBNInstance::random_device_());
