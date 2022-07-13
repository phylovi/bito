// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "gp_instance.hpp"

#include <chrono>
#include <iomanip>
#include <string>

#include "csv.hpp"
#include "driver.hpp"
#include "gp_operation.hpp"
#include "numerical_utils.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_probability.hpp"

using namespace GPOperations;  // NOLINT

void GPInstance::PrintStatus() {
  const auto tree_count = tree_collection_.TreeCount();
  const auto taxon_count = tree_collection_.TaxonCount();
  if (tree_count > 0) {
    std::cout << tree_count << " trees loaded on " << taxon_count << " leaves.\n";
  } else {
    std::cout << "No trees loaded.\n";
  }
  std::cout << alignment_.Data().size() << " sequences loaded.\n";
  std::cout << dag_.NodeCount() << " DAG nodes with "
            << dag_.EdgeCountWithLeafSubsplits() << " edges representing "
            << dag_.TopologyCount() << " trees.\n";
  std::cout << dag_.EdgeCountWithLeafSubsplits() << " continuous parameters.\n";
  if (HasEngine()) {
    std::cout << "Engine available using "
              << GetEngine()->GetPLVHandler().GetByteCount() / 1e9
              << "G virtual memory.\n";
  } else {
    std::cout << "Engine has not been made.\n";
  }
}

StringSizeMap GPInstance::DAGSummaryStatistics() { return dag_.SummaryStatistics(); }

void GPInstance::ReadFastaFile(const std::string &fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

void GPInstance::ReadNewickFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
}

void GPInstance::ReadNewickFileGZ(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFileGZ(fname));
}

void GPInstance::ReadNexusFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
}

void GPInstance::ReadNexusFileGZ(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFileGZ(fname));
}

void GPInstance::CheckSequencesLoaded() const {
  if (alignment_.SequenceCount() == 0) {
    Failwith(
        "Load an alignment into your GPInstance with which you wish to "
        "calculate phylogenetic likelihoods.");
  }
}

void GPInstance::CheckTreesLoaded() const {
  if (tree_collection_.TreeCount() == 0) {
    Failwith(
        "Load some trees into your GPInstance on which you wish to "
        "build your subsplit DAG.");
  }
}

void GPInstance::MakeDAG() {
  CheckTreesLoaded();
  dag_ = GPDAG(tree_collection_);
}

GPDAG &GPInstance::GetDAG() { return dag_; }

void GPInstance::PrintDAG() { dag_.Print(); }

void GPInstance::UseGradientOptimization(bool use_gradients) {
  use_gradients_ = use_gradients;
};

void GPInstance::MakeEngine(double rescaling_threshold) {
  CheckSequencesLoaded();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());

  MakeDAG();
  auto sbn_prior = dag_.BuildUniformOnTopologicalSupportPrior();
  auto unconditional_node_probabilities =
      dag_.UnconditionalNodeProbabilities(sbn_prior);
  auto inverted_sbn_prior =
      dag_.InvertedGPCSPProbabilities(sbn_prior, unconditional_node_probabilities);
  engine_ = std::make_unique<GPEngine>(
      std::move(site_pattern), dag_.NodeCountWithoutDAGRoot(), plv_count_per_node_,
      dag_.EdgeCountWithLeafSubsplits(), mmap_file_path_, rescaling_threshold,
      std::move(sbn_prior),
      unconditional_node_probabilities.segment(0, dag_.NodeCountWithoutDAGRoot()),
      std::move(inverted_sbn_prior), use_gradients_);
}

void GPInstance::ReinitializePriors() {
  auto sbn_prior = dag_.BuildUniformOnTopologicalSupportPrior();
  auto unconditional_node_probabilities =
      dag_.UnconditionalNodeProbabilities(sbn_prior);
  auto inverted_sbn_prior =
      dag_.InvertedGPCSPProbabilities(sbn_prior, unconditional_node_probabilities);
  engine_->InitializePriors(std::move(sbn_prior),
                            std::move(unconditional_node_probabilities.segment(
                                0, dag_.NodeCountWithoutDAGRoot())),
                            std::move(inverted_sbn_prior));
}

GPEngine *GPInstance::GetEngine() const {
  Assert(HasEngine(),
         "Engine not available. Call MakeEngine to make an engine for phylogenetic "
         "likelihood computation.");
  return engine_.get();
}

void GPInstance::ResizeEngineForDAG() {
  Assert(HasEngine(), "Engine not available. Call MakeEngine before resizing.");
  GetEngine()->GrowPLVs(GetDAG().NodeCountWithoutDAGRoot());
  GetEngine()->GrowGPCSPs(GetDAG().EdgeCountWithLeafSubsplits());
}

bool GPInstance::HasEngine() const { return engine_ != nullptr; }

void GPInstance::PrintEdgeIndexer() {
  std::cout << "Vector of taxon names: " << tree_collection_.TaxonNames() << std::endl;
  dag_.PrintEdgeIndexer();
}

void GPInstance::ProcessOperations(const GPOperationVector &operations) {
  GetEngine()->ProcessOperations(operations);
}

void GPInstance::ClearTreeCollectionAssociatedState() { dag_ = GPDAG(); }

void GPInstance::HotStartBranchLengths() {
  Assert(HasEngine(),
         "Please load and process some trees before calling HotStartBranchLengths.");
  GetEngine()->HotStartBranchLengths(tree_collection_, dag_.BuildEdgeIndexer());
}

SizeDoubleVectorMap GPInstance::GatherBranchLengths() {
  Assert(HasEngine(),
         "Please load and process some trees before calling GatherBranchLengths.");
  SizeDoubleVectorMap branch_lengths_from_sample =
      GetEngine()->GatherBranchLengths(tree_collection_, dag_.BuildEdgeIndexer());
  return branch_lengths_from_sample;
}

void GPInstance::TakeFirstBranchLength() {
  Assert(HasEngine(),
         "Please load and process some trees before calling TakeFirstBranchLength.");
  GetEngine()->TakeFirstBranchLength(tree_collection_, dag_.BuildEdgeIndexer());
}

void GPInstance::PopulatePLVs() { ProcessOperations(dag_.PopulatePLVs()); }

void GPInstance::ComputeLikelihoods() { ProcessOperations(dag_.ComputeLikelihoods()); }

void GPInstance::ComputeMarginalLikelihood() {
  ProcessOperations(dag_.MarginalLikelihood());
}

void GPInstance::EstimateBranchLengths(double tol, size_t max_iter, bool quiet,
                                       bool track_intermediate_iterations) {
  std::stringstream dev_null;

  auto &our_ostream = quiet ? dev_null : std::cout;
  auto now = std::chrono::high_resolution_clock::now;
  auto t_start = now();
  our_ostream << "Begin branch optimization\n";
  GPOperationVector branch_optimization_operations = dag_.BranchLengthOptimization();
  GPOperationVector marginal_lik_operations = dag_.MarginalLikelihood();
  GPOperationVector populate_plv_operations = dag_.PopulatePLVs();

  our_ostream << "Populating PLVs\n";
  PopulatePLVs();
  std::chrono::duration<double> warmup_duration = now() - t_start;
  t_start = now();
  our_ostream << "Computing initial likelihood\n";
  ProcessOperations(marginal_lik_operations);

  double current_marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();

  std::chrono::duration<double> initial_likelihood_duration = now() - t_start;
  t_start = now();

  for (size_t i = 0; i < max_iter; i++) {
    our_ostream << "Iteration: " << (i + 1) << std::endl;
    ProcessOperations(branch_optimization_operations);
    // #321 Replace with a cleaned up traversal.
    ProcessOperations(populate_plv_operations);
    ProcessOperations(marginal_lik_operations);

    double marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();

    if (track_intermediate_iterations) {
      our_ostream << "Tracking intermediate optimization values" << std::endl;
      IntermediateOptimizationValues();
    }

    our_ostream << "Current marginal log likelihood: ";
    our_ostream << std::setprecision(9) << current_marginal_log_lik << std::endl;
    our_ostream << "New marginal log likelihood: ";
    our_ostream << std::setprecision(9) << marginal_log_lik << std::endl;

    double avg_abs_change_perpcsp_branch_length =
        (GetEngine()->GetBranchLengthDifferences()).array().mean();

    our_ostream << "Average absolute change in branch lengths:";
    our_ostream << std::setprecision(9) << avg_abs_change_perpcsp_branch_length
                << std::endl;
    if (marginal_log_lik < current_marginal_log_lik) {
      our_ostream << "Marginal log likelihood decreased.\n";
    }
    if (avg_abs_change_perpcsp_branch_length < tol) {
      our_ostream << "Average absolute change in branch lengths converged. \n";
      break;
    }
    current_marginal_log_lik = marginal_log_lik;
  }

  std::chrono::duration<double> optimization_duration = now() - t_start;
  our_ostream << "\n# Timing Report\n";
  our_ostream << "warmup: " << warmup_duration.count() << "s\n";
  our_ostream << "initial likelihood: " << initial_likelihood_duration.count() << "s\n";
  our_ostream << "optimization: " << optimization_duration.count() << "s or "
              << optimization_duration.count() / 60 << "m\n";
}

void GPInstance::IntermediateOptimizationValues() {
  GPOperationVector compute_lik_operations = dag_.ComputeLikelihoods();
  ProcessOperations(compute_lik_operations);

  size_t col_idx = per_pcsp_branch_lengths_.cols() - 1;
  per_pcsp_branch_lengths_.col(col_idx) = GetEngine()->GetBranchLengths();
  per_pcsp_log_lik_.col(col_idx) = GetEngine()->GetPerGPCSPLogLikelihoods();

  per_pcsp_branch_lengths_.conservativeResize(Eigen::NoChange, col_idx + 2);
  per_pcsp_log_lik_.conservativeResize(Eigen::NoChange, col_idx + 2);
}

void GPInstance::EstimateSBNParameters() {
  std::cout << "Begin SBN parameter optimization\n";
  PopulatePLVs();
  ComputeLikelihoods();
  ProcessOperations(dag_.OptimizeSBNParameters());
}

void GPInstance::CalculateHybridMarginals() {
  std::cout << "Calculating hybrid marginals\n";
  PopulatePLVs();
  dag_.TopologicalEdgeTraversal([this](const NodeId parent_id,
                                       const bool is_edge_on_left,
                                       const NodeId child_id, const EdgeId edge_idx) {
    this->GetEngine()->ProcessQuartetHybridRequest(
        dag_.QuartetHybridRequestOf(parent_id, is_edge_on_left, child_id));
  });
}

size_t GPInstance::GetEdgeIndexForLeafNode(const Bitset &parent_subsplit,
                                           const Node *leaf_node) const {
  Assert(leaf_node->IsLeaf(), "Only leaf node is permitted.");
  return dag_.GetEdgeIdx(parent_subsplit,
                         Bitset::LeafSubsplitOfNonemptyClade(leaf_node->Leaves()));
}

RootedTreeCollection GPInstance::TreesWithGPBranchLengthsOfTopologies(
    Node::NodePtrVec &&topologies) const {
  const EigenVectorXd gpcsp_indexed_branch_lengths = GetEngine()->GetBranchLengths();
  RootedTree::RootedTreeVector tree_vector;

  for (auto &root_node : topologies) {
    size_t node_count = 2 * root_node->LeafCount() - 1;
    std::vector<double> branch_lengths(node_count);

    root_node->RootedPCSPPreorder(
        [this, &branch_lengths, &gpcsp_indexed_branch_lengths](
            const Node *sister, const Node *focal, const Node *child0,
            const Node *child1) {
          Bitset parent_subsplit = Bitset::Subsplit(sister->Leaves(), focal->Leaves());
          Bitset child_subsplit = Bitset::Subsplit(child0->Leaves(), child1->Leaves());
          size_t gpcsp_idx = dag_.GetEdgeIdx(parent_subsplit, child_subsplit);
          branch_lengths[focal->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];

          if (sister->IsLeaf()) {
            gpcsp_idx = GetEdgeIndexForLeafNode(parent_subsplit, sister);
            branch_lengths[sister->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];
          }
          if (child0->IsLeaf()) {
            gpcsp_idx = GetEdgeIndexForLeafNode(child_subsplit, child0);
            branch_lengths[child0->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];
          }
          if (child1->IsLeaf()) {
            gpcsp_idx = GetEdgeIndexForLeafNode(child_subsplit, child1);
            branch_lengths[child1->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];
          }
        },
        false);

    tree_vector.emplace_back(root_node, std::move(branch_lengths));
  }

  return RootedTreeCollection(tree_vector, tree_collection_.TagTaxonMap());
}

RootedTreeCollection GPInstance::GenerateCompleteRootedTreeCollection() {
  return TreesWithGPBranchLengthsOfTopologies(dag_.GenerateAllTopologies());
}

void GPInstance::GetPerGPCSPLogLikelihoodSurfaces(size_t steps, double scale_min,
                                                  double scale_max) {
  const EigenVectorXd optimized_branch_lengths = GetEngine()->GetBranchLengths();

  size_t gpcsp_count = optimized_branch_lengths.size();
  const EigenVectorXd scaling_vector =
      EigenVectorXd::LinSpaced(steps, scale_min, scale_max);
  per_pcsp_lik_surfaces_ = EigenMatrixXd(gpcsp_count * steps, 2);

  for (EdgeId gpcsp_idx = 0; gpcsp_idx < gpcsp_count; gpcsp_idx++) {
    EigenVectorXd gpcsp_new_branch_lengths =
        scaling_vector * optimized_branch_lengths[gpcsp_idx];
    EigenVectorXd new_branch_length_vector = optimized_branch_lengths;

    for (size_t i = 0; i < steps; i++) {
      new_branch_length_vector[gpcsp_idx] = gpcsp_new_branch_lengths[i];
      GetEngine()->SetBranchLengths(new_branch_length_vector);
      PopulatePLVs();
      ComputeLikelihoods();

      size_t matrix_position = gpcsp_count * i;
      per_pcsp_lik_surfaces_(matrix_position + gpcsp_idx, 0) =
          gpcsp_new_branch_lengths[i];
      per_pcsp_lik_surfaces_(matrix_position + gpcsp_idx, 1) =
          GetEngine()->GetPerGPCSPLogLikelihoods(gpcsp_idx, 1)(0, 0);
    }
  }
  // Reset back to optimized branch lengths
  GetEngine()->SetBranchLengths(optimized_branch_lengths);
}

void GPInstance::PerturbAndTrackValuesFromOptimization() {
  const EigenVectorXd optimized_branch_lengths = GetEngine()->GetBranchLengths();
  const EigenVectorXd optimized_per_pcsp_llhs =
      GetEngine()->GetPerGPCSPLogLikelihoods();

  size_t gpcsp_count = optimized_branch_lengths.size();
  EigenVectorXi run_counts = EigenVectorXi::Zero(gpcsp_count);
  EigenMatrixXd tracked_optimization_values(1, 2);

  const auto pretty_indexer = PrettyIndexer();
  StringVector pretty_index_vector;
  GPOperationVector branch_optimization_operations = dag_.BranchLengthOptimization();

  for (size_t gpcsp_idx = 0; gpcsp_idx < gpcsp_count; gpcsp_idx++) {
    double optimized_llh = optimized_per_pcsp_llhs[gpcsp_idx];
    double current_branch_length = 0.1;

    while (true) {
      run_counts[gpcsp_idx]++;

      EigenVectorXd new_branch_length_vector = optimized_branch_lengths;
      new_branch_length_vector[gpcsp_idx] = current_branch_length;
      GetEngine()->SetBranchLengths(new_branch_length_vector);

      PopulatePLVs();
      ComputeLikelihoods();

      double current_llh = GetEngine()->GetPerGPCSPLogLikelihoods(gpcsp_idx, 1)(0, 0);

      tracked_optimization_values.row(tracked_optimization_values.rows() - 1)
          << current_branch_length,
          current_llh;
      tracked_optimization_values.conservativeResize(
          tracked_optimization_values.rows() + 1, Eigen::NoChange);

      if (fabs(current_llh - optimized_llh) < 1e-3 || run_counts[gpcsp_idx] > 5) {
        break;
      } else {
        ProcessOperations(branch_optimization_operations);
        current_branch_length = GetEngine()->GetBranchLengths()[gpcsp_idx];
      }
    }
    pretty_index_vector.insert(pretty_index_vector.end(), run_counts[gpcsp_idx],
                               pretty_indexer.at(gpcsp_idx));
  }
  // Reset back to optimized branch lengths
  GetEngine()->SetBranchLengths(optimized_branch_lengths);

  tracked_optimization_values.conservativeResize(tracked_optimization_values.rows() - 1,
                                                 Eigen::NoChange);

  tracked_values_after_perturbing_.reserve(tracked_optimization_values.rows());
  for (int i = 0; i < tracked_optimization_values.rows(); i++) {
    tracked_values_after_perturbing_.push_back(
        {pretty_index_vector.at(i), tracked_optimization_values.row(i)});
  }
}

StringVector GPInstance::PrettyIndexer() const {
  StringVector pretty_representation(dag_.BuildEdgeIndexer().size());
  // #350 consider use of edge vs pcsp here.
  for (const auto &[edge, idx] : dag_.BuildEdgeIndexer()) {
    pretty_representation[idx] = edge.PCSPToString();
  }
  return pretty_representation;
}

StringDoubleVector GPInstance::PrettyIndexedVector(EigenConstVectorXdRef v) {
  StringDoubleVector result;
  result.reserve(v.size());
  const auto pretty_indexer = PrettyIndexer();
  Assert(v.size() <= static_cast<Eigen::Index>(pretty_indexer.size()),
         "v is too long in PrettyIndexedVector");
  for (Eigen::Index i = 0; i < v.size(); i++) {
    result.push_back({pretty_indexer.at(i), v(i)});
  }
  return result;
}

VectorOfStringAndEigenVectorXdPairs GPInstance::PrettyIndexedMatrix(
    EigenConstMatrixXdRef m) {
  VectorOfStringAndEigenVectorXdPairs result;
  result.reserve(m.rows());
  const auto pretty_indexer = PrettyIndexer();
  for (int i = 0; i < m.rows(); i++) {
    int idx = i % pretty_indexer.size();
    result.push_back({pretty_indexer.at(idx), m.row(i)});
  }
  return result;
}

EigenConstVectorXdRef GPInstance::GetSBNParameters() {
  return GetEngine()->GetSBNParameters();
}

StringDoubleVector GPInstance::PrettyIndexedSBNParameters() {
  return PrettyIndexedVector(GetSBNParameters());
}

StringDoubleVector GPInstance::PrettyIndexedBranchLengths() {
  return PrettyIndexedVector(GetEngine()->GetBranchLengths());
}

StringDoubleVector GPInstance::PrettyIndexedPerGPCSPLogLikelihoods() {
  return PrettyIndexedVector(GetEngine()->GetPerGPCSPLogLikelihoods());
}

StringDoubleVector GPInstance::PrettyIndexedPerGPCSPComponentsOfFullLogMarginal() {
  return PrettyIndexedVector(GetEngine()->GetPerGPCSPComponentsOfFullLogMarginal());
}

VectorOfStringAndEigenVectorXdPairs
GPInstance::PrettyIndexedIntermediateBranchLengths() {
  return PrettyIndexedMatrix(per_pcsp_branch_lengths_);
}

VectorOfStringAndEigenVectorXdPairs
GPInstance::PrettyIndexedIntermediatePerGPCSPLogLikelihoods() {
  return PrettyIndexedMatrix(per_pcsp_log_lik_);
}

VectorOfStringAndEigenVectorXdPairs
GPInstance::PrettyIndexedPerGPCSPLogLikelihoodSurfaces() {
  return PrettyIndexedMatrix(per_pcsp_lik_surfaces_);
}

void GPInstance::SBNParametersToCSV(const std::string &file_path) {
  CSV::StringDoubleVectorToCSV(PrettyIndexedSBNParameters(), file_path);
}

void GPInstance::SBNPriorToCSV(const std::string &file_path) {
  CSV::StringDoubleVectorToCSV(
      PrettyIndexedVector(dag_.BuildUniformOnTopologicalSupportPrior()), file_path);
}

void GPInstance::BranchLengthsToCSV(const std::string &file_path) {
  CSV::StringDoubleVectorToCSV(PrettyIndexedBranchLengths(), file_path);
}

void GPInstance::PerGPCSPLogLikelihoodsToCSV(const std::string &file_path) {
  CSV::StringDoubleVectorToCSV(PrettyIndexedPerGPCSPLogLikelihoods(), file_path);
}

void GPInstance::PerPCSPIndexedMatrixToCSV(
    VectorOfStringAndEigenVectorXdPairs per_pcsp_indexed_matrix,
    const std::string &file_path) {
  std::ofstream out_stream(file_path);

  for (const auto &[s, eigen] : per_pcsp_indexed_matrix) {
    out_stream << s;
    for (const auto &value : eigen) {
      out_stream << "," << std::setprecision(9) << value;
    }
    out_stream << std::endl;
  }
  if (out_stream.bad()) {
    Failwith("Failure writing to " + file_path);
  }
  out_stream.close();
}

void GPInstance::IntermediateBranchLengthsToCSV(const std::string &file_path) {
  return PerPCSPIndexedMatrixToCSV(PrettyIndexedIntermediateBranchLengths(), file_path);
}

void GPInstance::IntermediatePerGPCSPLogLikelihoodsToCSV(const std::string &file_path) {
  return PerPCSPIndexedMatrixToCSV(PrettyIndexedIntermediatePerGPCSPLogLikelihoods(),
                                   file_path);
}

void GPInstance::PerGPCSPLogLikelihoodSurfacesToCSV(const std::string &file_path) {
  std::ofstream out_stream(file_path);
  VectorOfStringAndEigenVectorXdPairs vect =
      PrettyIndexedPerGPCSPLogLikelihoodSurfaces();

  for (const auto &[s, eigen] : vect) {
    out_stream << s;
    for (const auto &value : eigen) {
      out_stream << "," << std::setprecision(9) << value;
    }
    out_stream << std::endl;
  }
  if (out_stream.bad()) {
    Failwith("Failure writing to " + file_path);
  }
  out_stream.close();
}

void GPInstance::TrackedOptimizationValuesToCSV(const std::string &file_path) {
  return PerPCSPIndexedMatrixToCSV(tracked_values_after_perturbing_, file_path);
}

RootedTreeCollection GPInstance::CurrentlyLoadedTreesWithGPBranchLengths() {
  Node::NodePtrVec topologies;
  for (const auto &tree : tree_collection_.Trees()) {
    topologies.push_back(tree.Topology()->DeepCopy());
  }
  return TreesWithGPBranchLengthsOfTopologies(std::move(topologies));
}

RootedTreeCollection GPInstance::CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths(
    const std::string &pcsp_string) {
  const BitsetSizeMap &indexer = dag_.BuildEdgeIndexer();
  Bitset pcsp(pcsp_string);
  auto search = indexer.find(pcsp);
  if (search == indexer.end()) {
    Failwith("Don't have " + pcsp_string + " as a PCSP in the instance!");
  }
  auto pcsp_index = search->second;

  Node::NodePtrVec topologies;
  for (const auto &tree : tree_collection_.Trees()) {
    auto indexer_representation = dag_.IndexerRepresentationOf(
        indexer, tree.Topology(), std::numeric_limits<size_t>::max());
    if (std::find(indexer_representation.begin(), indexer_representation.end(),
                  pcsp_index) != indexer_representation.end()) {
      topologies.push_back(tree.Topology()->DeepCopy());
    }
  }
  return TreesWithGPBranchLengthsOfTopologies(std::move(topologies));
}

void GPInstance::ExportTrees(const std::string &out_path) {
  auto trees = CurrentlyLoadedTreesWithGPBranchLengths();
  trees.ToNewickFile(out_path);
}

void GPInstance::ExportTreesWithAPCSP(const std::string &pcsp_string,
                                      const std::string &out_path) {
  auto trees = CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths(pcsp_string);
  trees.ToNewickFile(out_path);
}

void GPInstance::ExportAllGeneratedTrees(const std::string &out_path) {
  auto trees = GenerateCompleteRootedTreeCollection();
  trees.ToNewickFile(out_path);
}

void GPInstance::ExportAllGeneratedTopologies(const std::string &out_path) {
  TreeCollection::UnitBranchLengthTreesOf(dag_.GenerateAllTopologies(),
                                          tree_collection_.TagTaxonMap())
      .ToNewickTopologyFile(out_path);
}

void GPInstance::LoadAllGeneratedTrees() {
  tree_collection_ = GenerateCompleteRootedTreeCollection();
}

StringVector GPInstance::GetTaxonNames() const { return tree_collection_.TaxonNames(); }

EigenVectorXd GPInstance::GetBranchLengths() const {
  return GetEngine()->GetBranchLengths();
}

void GPInstance::SubsplitDAGToDot(const std::string &out_path,
                                  bool show_index_labels) const {
  std::ofstream out_stream(out_path);
  out_stream << dag_.ToDot(show_index_labels) << std::endl;
  if (out_stream.bad()) {
    Failwith("Failure writing to " + out_path);
  }
  out_stream.close();
}

void GPInstance::MakeNNIEngine() {
  nni_engine_ = std::make_unique<NNIEngine>(dag_, *engine_);
}

NNIEngine &GPInstance::GetNNIEngine() {
  Assert(nni_engine_, "GPInstance::GetNNIEngine() when nni_engine has not been made.");
  return *nni_engine_;
}
