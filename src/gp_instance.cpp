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
#include "stopwatch.hpp"

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
  std::cout << GetDAG().NodeCount() << " DAG nodes with "
            << GetDAG().EdgeCountWithLeafSubsplits() << " edges representing "
            << GetDAG().TopologyCount() << " trees.\n";
  std::cout << GetDAG().EdgeCountWithLeafSubsplits() << " continuous parameters.\n";
  if (HasGPEngine()) {
    std::cout << "Engine available using "
              << GetGPEngine().GetPLVHandler().GetByteCount() / 1e9
              << "G virtual memory.\n";
  } else {
    std::cout << "Engine has not been made.\n";
  }
}

StringSizeMap GPInstance::DAGSummaryStatistics() {
  return GetDAG().SummaryStatistics();
}

void GPInstance::ReadFastaFile(const std::string &fname) {
  alignment_ = Alignment::ReadFasta(fname);
  fasta_path_ = fname;
}

void GPInstance::ReadNewickFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
  newick_path_ = fname;
}

void GPInstance::ReadNewickFileGZ(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFileGZ(fname));
  newick_path_ = fname;
}

void GPInstance::ReadNexusFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
  nexus_path_ = fname;
}

void GPInstance::ReadNexusFileGZ(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFileGZ(fname));
  nexus_path_ = fname;
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
  dag_ = std::make_unique<GPDAG>(tree_collection_);
}

GPDAG &GPInstance::GetDAG() {
  Assert(HasDAG(), "DAG not available. Call MakeDAG.");
  return *dag_.get();
}

const GPDAG &GPInstance::GetDAG() const {
  Assert(HasDAG(), "DAG not available. Call MakeDAG.");
  return *dag_.get();
}

bool GPInstance::HasDAG() const { return dag_ != nullptr; }

void GPInstance::PrintDAG() { GetDAG().Print(); }

SitePattern GPInstance::MakeSitePattern() const {
  CheckSequencesLoaded();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
  return site_pattern;
}

// ** GP Engine

void GPInstance::MakeGPEngine(double rescaling_threshold, bool use_gradients) {
  std::string mmap_gp_path = mmap_file_path_.value() + ".gp";
  auto site_pattern = MakeSitePattern();
  if (!HasDAG()) {
    MakeDAG();
  }
  auto sbn_prior = GetDAG().BuildUniformOnTopologicalSupportPrior();
  auto unconditional_node_probabilities =
      GetDAG().UnconditionalNodeProbabilities(sbn_prior);
  auto inverted_sbn_prior =
      GetDAG().InvertedGPCSPProbabilities(sbn_prior, unconditional_node_probabilities);

  gp_engine_ = std::make_unique<GPEngine>(
      std::move(site_pattern), GetDAG().NodeCountWithoutDAGRoot(),
      GetDAG().EdgeCountWithLeafSubsplits(), mmap_gp_path, rescaling_threshold,
      std::move(sbn_prior),
      unconditional_node_probabilities.segment(0, GetDAG().NodeCountWithoutDAGRoot()),
      std::move(inverted_sbn_prior), use_gradients);
}

void GPInstance::ReinitializePriors() {
  auto sbn_prior = GetDAG().BuildUniformOnTopologicalSupportPrior();
  auto unconditional_node_probabilities =
      GetDAG().UnconditionalNodeProbabilities(sbn_prior);
  auto inverted_sbn_prior =
      GetDAG().InvertedGPCSPProbabilities(sbn_prior, unconditional_node_probabilities);
  GetGPEngine().InitializePriors(std::move(sbn_prior),
                                 std::move(unconditional_node_probabilities.segment(
                                     0, GetDAG().NodeCountWithoutDAGRoot())),
                                 std::move(inverted_sbn_prior));
}

GPEngine &GPInstance::GetGPEngine() const {
  Assert(HasGPEngine(),
         "Engine not available. Call MakeGPEngine to make an engine for phylogenetic "
         "likelihood computation.");
  return *gp_engine_.get();
}

void GPInstance::ResizeEngineForDAG() {
  Assert(HasGPEngine(), "Engine not available. Call MakeGPEngine before resizing.");
  GetGPEngine().GrowPLVs(GetDAG().NodeCountWithoutDAGRoot());
  GetGPEngine().GrowGPCSPs(GetDAG().EdgeCountWithLeafSubsplits());
}

bool GPInstance::HasGPEngine() const { return gp_engine_ != nullptr; }

void GPInstance::PrintEdgeIndexer() {
  std::cout << "Vector of taxon names: " << tree_collection_.TaxonNames() << std::endl;
  GetDAG().PrintEdgeIndexer();
}

void GPInstance::ProcessOperations(const GPOperationVector &operations) {
  GetGPEngine().ProcessOperations(operations);
}

void GPInstance::ClearTreeCollectionAssociatedState() { GetDAG() = GPDAG(); }

void GPInstance::HotStartBranchLengths() {
  Assert(HasGPEngine(),
         "Please load and process some trees before calling HotStartBranchLengths.");
  GetGPEngine().HotStartBranchLengths(tree_collection_, GetDAG().BuildEdgeIndexer());
}

SizeDoubleVectorMap GPInstance::GatherBranchLengths() {
  Assert(HasGPEngine(),
         "Please load and process some trees before calling GatherBranchLengths.");
  SizeDoubleVectorMap branch_lengths_from_sample =
      GetGPEngine().GatherBranchLengths(tree_collection_, GetDAG().BuildEdgeIndexer());
  return branch_lengths_from_sample;
}

void GPInstance::TakeFirstBranchLength() {
  Assert(HasGPEngine(),
         "Please load and process some trees before calling TakeFirstBranchLength.");
  GetGPEngine().TakeFirstBranchLength(tree_collection_, GetDAG().BuildEdgeIndexer());
}

void GPInstance::PopulatePLVs() { ProcessOperations(GetDAG().PopulatePLVs()); }

void GPInstance::ComputeLikelihoods() {
  ProcessOperations(GetDAG().ComputeLikelihoods());
}

void GPInstance::ComputeMarginalLikelihood() {
  ProcessOperations(GetDAG().MarginalLikelihood());
}

void GPInstance::EstimateBranchLengths(double tol, size_t max_iter, bool quiet,
                                       bool track_intermediate_iterations,
                                       std::optional<OptimizationMethod> method) {
  std::stringstream dev_null;
  auto &our_ostream = quiet ? dev_null : std::cout;
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);

  if (method.has_value()) {
    GetGPEngine().SetOptimizationMethod(method.value());
  }
  GetGPEngine().ResetOptimizationCount();

  our_ostream << "Begin branch optimization\n";
  // GetDAG().ReinitializeTidyVectors();
  GPOperationVector branch_optimization_operations =
      GetDAG().BranchLengthOptimization();
  GPOperationVector marginal_lik_operations = GetDAG().MarginalLikelihood();
  GPOperationVector populate_plv_operations = GetDAG().PopulatePLVs();

  our_ostream << "Populating PLVs\n";
  PopulatePLVs();
  auto warmup_duration = timer.Lap();
  our_ostream << "Computing initial likelihood\n";
  ProcessOperations(marginal_lik_operations);
  double current_marginal_log_lik = GetGPEngine().GetLogMarginalLikelihood();
  auto initial_likelihood_duration = timer.Lap();

  for (size_t i = 0; i < max_iter; i++) {
    our_ostream << "Iteration: " << (i + 1) << std::endl;
    ProcessOperations(branch_optimization_operations);
    // #321 Replace with a cleaned up traversal.
    ProcessOperations(populate_plv_operations);
    ProcessOperations(marginal_lik_operations);
    double marginal_log_lik = GetGPEngine().GetLogMarginalLikelihood();

    if (track_intermediate_iterations) {
      our_ostream << "Tracking intermediate optimization values" << std::endl;
      IntermediateOptimizationValues();
    }

    our_ostream << "Current marginal log likelihood: ";
    our_ostream << std::setprecision(9) << current_marginal_log_lik << std::endl;
    our_ostream << "New marginal log likelihood: ";
    our_ostream << std::setprecision(9) << marginal_log_lik << std::endl;
    double avg_abs_change_perpcsp_branch_length =
        GetGPEngine().GetBranchLengthDifferences().array().mean();
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
    GetGPEngine().IncrementOptimizationCount();
  }

  timer.Stop();
  auto optimization_duration = timer.GetTotal();
  our_ostream << "\n# Timing Report\n";
  our_ostream << "warmup: " << warmup_duration << "s\n";
  our_ostream << "initial likelihood: " << initial_likelihood_duration << "s\n";
  our_ostream << "optimization: " << optimization_duration << "s or "
              << optimization_duration / 60 << "m\n";
}

void GPInstance::EstimateTPBranchLengths(double tol, size_t max_iter, bool quiet,
                                         bool track_intermediate_iterations,
                                         std::optional<OptimizationMethod> method) {
  std::stringstream dev_null;
  auto &our_ostream = quiet ? dev_null : std::cout;

  auto &tp_engine = GetTPEngine();
  auto &dag = GetDAG();
  auto &branch_handler = GetTPEngine().GetLikelihoodEvalEngine().GetDAGBranchHandler();
  auto &diffs = branch_handler.GetBranchDifferences().GetData();
  auto &branches = branch_handler.GetBranchLengths().GetData();
  if (method.has_value()) {
    tp_engine.GetDAGBranchHandler().SetOptimizationMethod(method.value());
  }
  tp_engine.GetLikelihoodEvalEngine().ResetOptimizationCount();

  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);

  our_ostream << "Begin branch optimization\n";
  tp_engine.GetLikelihoodEvalEngine().InitializeBranchLengthHandler();

  our_ostream << "Computing Likelihoods\n";
  tp_engine.GetLikelihoodEvalEngine().Initialize();
  tp_engine.GetLikelihoodEvalEngine().ComputeScores();
  our_ostream << "Likelihoods: " << std::endl
              << tp_engine.GetTopTreeLikelihoods() << std::endl;
  auto cur_likelihoods = branches.segment(1, dag.EdgeCountWithLeafSubsplits()).mean();
  auto prv_likelihoods = cur_likelihoods;
  auto warmup_duration = timer.Lap();
  auto initial_likelihood_duration = timer.Lap();

  for (size_t i = 0; i < max_iter; i++) {
    tp_engine.GetLikelihoodEvalEngine().BranchLengthOptimization();
    tp_engine.GetLikelihoodEvalEngine().ComputeScores();
    auto cur_likelihoods = tp_engine.GetTopTreeLikelihoods()
                               .segment(0, GetDAG().EdgeCountWithLeafSubsplits())
                               .mean();

    our_ostream << "Iteration: " << (i + 1) << std::endl;

    our_ostream << "Previous likelihoods: ";
    our_ostream << std::setprecision(9) << prv_likelihoods << std::endl;
    our_ostream << "Current likelihoods: ";
    our_ostream << std::setprecision(9) << cur_likelihoods << std::endl;

    double avg_abs_change_perpcsp_branch_length =
        diffs.segment(1, dag.EdgeCountWithLeafSubsplits()).mean();

    our_ostream << "Average absolute change in branch lengths:";
    our_ostream << std::setprecision(9) << avg_abs_change_perpcsp_branch_length
                << std::endl;
    if (prv_likelihoods < cur_likelihoods) {
      our_ostream << "Marginal log likelihood decreased.\n";
    }
    if (avg_abs_change_perpcsp_branch_length < tol) {
      our_ostream << "Average absolute change in branch lengths converged. \n";
      break;
    }
    prv_likelihoods = cur_likelihoods;

    tp_engine.GetLikelihoodEvalEngine().IncrementOptimizationCount();
  }

  auto optimization_duration = timer.GetTotal();
  our_ostream << "\n# Timing Report\n";
  our_ostream << "warmup: " << warmup_duration << "s\n";
  our_ostream << "initial likelihood: " << initial_likelihood_duration << "s\n";
  our_ostream << "optimization: " << optimization_duration << "s or "
              << optimization_duration / 60 << "m\n";
}

void GPInstance::SetOptimizationMethod(const OptimizationMethod method) {
  GetGPEngine().SetOptimizationMethod(method);
}

void GPInstance::UseGradientOptimization(const bool use_gradients) {
  GetGPEngine().UseGradientOptimization(use_gradients);
}

void GPInstance::IntermediateOptimizationValues() {
  GPOperationVector compute_lik_operations = GetDAG().ComputeLikelihoods();
  ProcessOperations(compute_lik_operations);

  size_t col_idx = per_pcsp_branch_lengths_.cols() - 1;
  per_pcsp_branch_lengths_.col(col_idx) = GetGPEngine().GetBranchLengths();
  per_pcsp_log_lik_.col(col_idx) = GetGPEngine().GetPerGPCSPLogLikelihoods();

  per_pcsp_branch_lengths_.conservativeResize(Eigen::NoChange, col_idx + 2);
  per_pcsp_log_lik_.conservativeResize(Eigen::NoChange, col_idx + 2);
}

void GPInstance::EstimateSBNParameters() {
  std::cout << "Begin SBN parameter optimization\n";
  PopulatePLVs();
  ComputeLikelihoods();
  ProcessOperations(GetDAG().OptimizeSBNParameters());
}

void GPInstance::CalculateHybridMarginals() {
  std::cout << "Calculating hybrid marginals\n";
  PopulatePLVs();
  GetDAG().TopologicalEdgeTraversal(
      [this](const NodeId parent_id, const bool is_edge_on_left, const NodeId child_id,
             const EdgeId edge_idx) {
        this->GetGPEngine().ProcessQuartetHybridRequest(
            GetDAG().QuartetHybridRequestOf(parent_id, is_edge_on_left, child_id));
      });
}

EdgeId GPInstance::GetEdgeIndexForLeafNode(const Bitset &parent_subsplit,
                                           const Node *leaf_node) const {
  Assert(leaf_node->IsLeaf(), "Only leaf node is permitted.");
  return GetDAG().GetEdgeIdx(parent_subsplit,
                             Bitset::LeafSubsplitOfNonemptyClade(leaf_node->Leaves()));
}

RootedTreeCollection GPInstance::TreesWithGPBranchLengthsOfTopologies(
    Node::NodePtrVec &&topologies) const {
  const EigenVectorXd gpcsp_indexed_branch_lengths = GetGPEngine().GetBranchLengths();
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
          EdgeId gpcsp_idx = GetDAG().GetEdgeIdx(parent_subsplit, child_subsplit);
          branch_lengths[focal->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx.value_];

          if (sister->IsLeaf()) {
            gpcsp_idx = GetEdgeIndexForLeafNode(parent_subsplit, sister);
            branch_lengths[sister->Id()] =
                gpcsp_indexed_branch_lengths[gpcsp_idx.value_];
          }
          if (child0->IsLeaf()) {
            gpcsp_idx = GetEdgeIndexForLeafNode(child_subsplit, child0);
            branch_lengths[child0->Id()] =
                gpcsp_indexed_branch_lengths[gpcsp_idx.value_];
          }
          if (child1->IsLeaf()) {
            gpcsp_idx = GetEdgeIndexForLeafNode(child_subsplit, child1);
            branch_lengths[child1->Id()] =
                gpcsp_indexed_branch_lengths[gpcsp_idx.value_];
          }
        },
        false);

    tree_vector.emplace_back(root_node, std::move(branch_lengths));
  }

  return RootedTreeCollection(tree_vector, tree_collection_.TagTaxonMap());
}

RootedTreeCollection GPInstance::GenerateCompleteRootedTreeCollection() {
  return TreesWithGPBranchLengthsOfTopologies(GetDAG().GenerateAllTopologies());
}

void GPInstance::GetPerGPCSPLogLikelihoodSurfaces(size_t steps, double scale_min,
                                                  double scale_max) {
  const EigenVectorXd optimized_branch_lengths = GetGPEngine().GetBranchLengths();

  size_t gpcsp_count = optimized_branch_lengths.size();
  const EigenVectorXd scaling_vector =
      EigenVectorXd::LinSpaced(steps, scale_min, scale_max);
  per_pcsp_lik_surfaces_ = EigenMatrixXd(gpcsp_count * steps, 2);

  for (EdgeId gpcsp_idx = EdgeId(0); gpcsp_idx < gpcsp_count; gpcsp_idx++) {
    EigenVectorXd gpcsp_new_branch_lengths =
        scaling_vector * optimized_branch_lengths[gpcsp_idx.value_];
    EigenVectorXd new_branch_length_vector = optimized_branch_lengths;

    for (size_t i = 0; i < steps; i++) {
      new_branch_length_vector[gpcsp_idx.value_] = gpcsp_new_branch_lengths[i];
      GetGPEngine().SetBranchLengths(new_branch_length_vector);
      PopulatePLVs();
      ComputeLikelihoods();

      size_t matrix_position = gpcsp_count * i;
      per_pcsp_lik_surfaces_(matrix_position + gpcsp_idx.value_, 0) =
          gpcsp_new_branch_lengths[i];
      per_pcsp_lik_surfaces_(matrix_position + gpcsp_idx.value_, 1) =
          GetGPEngine().GetPerGPCSPLogLikelihoods(gpcsp_idx.value_, 1)(0, 0);
    }
  }
  // Reset back to optimized branch lengths
  GetGPEngine().SetBranchLengths(optimized_branch_lengths);
}

void GPInstance::PerturbAndTrackValuesFromOptimization() {
  const EigenVectorXd optimized_branch_lengths = GetGPEngine().GetBranchLengths();
  const EigenVectorXd optimized_per_pcsp_llhs =
      GetGPEngine().GetPerGPCSPLogLikelihoods();

  size_t gpcsp_count = optimized_branch_lengths.size();
  EigenVectorXi run_counts = EigenVectorXi::Zero(gpcsp_count);
  EigenMatrixXd tracked_optimization_values(1, 2);

  const auto pretty_indexer = PrettyIndexer();
  StringVector pretty_index_vector;
  GPOperationVector branch_optimization_operations =
      GetDAG().BranchLengthOptimization();

  for (size_t gpcsp_idx = 0; gpcsp_idx < gpcsp_count; gpcsp_idx++) {
    double optimized_llh = optimized_per_pcsp_llhs[gpcsp_idx];
    double current_branch_length = 0.1;

    while (true) {
      run_counts[gpcsp_idx]++;

      EigenVectorXd new_branch_length_vector = optimized_branch_lengths;
      new_branch_length_vector[gpcsp_idx] = current_branch_length;
      GetGPEngine().SetBranchLengths(new_branch_length_vector);

      PopulatePLVs();
      ComputeLikelihoods();

      double current_llh = GetGPEngine().GetPerGPCSPLogLikelihoods(gpcsp_idx, 1)(0, 0);

      tracked_optimization_values.row(tracked_optimization_values.rows() - 1)
          << current_branch_length,
          current_llh;
      tracked_optimization_values.conservativeResize(
          tracked_optimization_values.rows() + 1, Eigen::NoChange);

      if (fabs(current_llh - optimized_llh) < 1e-3 || run_counts[gpcsp_idx] > 5) {
        break;
      } else {
        ProcessOperations(branch_optimization_operations);
        current_branch_length = GetGPEngine().GetBranchLengths()[gpcsp_idx];
      }
    }
    pretty_index_vector.insert(pretty_index_vector.end(), run_counts[gpcsp_idx],
                               pretty_indexer.at(gpcsp_idx));
  }
  // Reset back to optimized branch lengths
  GetGPEngine().SetBranchLengths(optimized_branch_lengths);

  tracked_optimization_values.conservativeResize(tracked_optimization_values.rows() - 1,
                                                 Eigen::NoChange);

  tracked_values_after_perturbing_.reserve(tracked_optimization_values.rows());
  for (int i = 0; i < tracked_optimization_values.rows(); i++) {
    tracked_values_after_perturbing_.push_back(
        {pretty_index_vector.at(i), tracked_optimization_values.row(i)});
  }
}

StringVector GPInstance::PrettyIndexer() const {
  StringVector pretty_representation(GetDAG().BuildEdgeIndexer().size());
  // #350 consider use of edge vs pcsp here.
  for (const auto &[edge, idx] : GetDAG().BuildEdgeIndexer()) {
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
  return GetGPEngine().GetSBNParameters();
}

StringDoubleVector GPInstance::PrettyIndexedSBNParameters() {
  return PrettyIndexedVector(GetSBNParameters());
}

StringDoubleVector GPInstance::PrettyIndexedBranchLengths() {
  return PrettyIndexedVector(GetGPEngine().GetBranchLengths());
}

StringDoubleVector GPInstance::PrettyIndexedPerGPCSPLogLikelihoods() {
  return PrettyIndexedVector(GetGPEngine().GetPerGPCSPLogLikelihoods());
}

StringDoubleVector GPInstance::PrettyIndexedPerGPCSPComponentsOfFullLogMarginal() {
  return PrettyIndexedVector(GetGPEngine().GetPerGPCSPComponentsOfFullLogMarginal());
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
      PrettyIndexedVector(GetDAG().BuildUniformOnTopologicalSupportPrior()), file_path);
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
  const BitsetSizeMap &indexer = GetDAG().BuildEdgeIndexer();
  Bitset pcsp(pcsp_string);
  auto search = indexer.find(pcsp);
  if (search == indexer.end()) {
    Failwith("Don't have " + pcsp_string + " as a PCSP in the instance!");
  }
  auto pcsp_index = search->second;

  Node::NodePtrVec topologies;
  for (const auto &tree : tree_collection_.Trees()) {
    auto indexer_representation = GetDAG().IndexerRepresentationOf(
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
  TreeCollection::UnitBranchLengthTreesOf(GetDAG().GenerateAllTopologies(),
                                          tree_collection_.TagTaxonMap())
      .ToNewickTopologyFile(out_path);
}

void GPInstance::LoadAllGeneratedTrees() {
  tree_collection_ = GenerateCompleteRootedTreeCollection();
}

StringVector GPInstance::GetTaxonNames() const { return tree_collection_.TaxonNames(); }

EigenVectorXd GPInstance::GetBranchLengths() const {
  return GetGPEngine().GetBranchLengths();
}

EigenVectorXd GPInstance::GetPerPCSPLogLikelihoods() const {
  return GetGPEngine().GetPerGPCSPLogLikelihoods();
}

void GPInstance::SubsplitDAGToDot(const std::string &out_path,
                                  bool show_index_labels) const {
  std::ofstream out_stream(out_path);
  out_stream << GetDAG().ToDot(show_index_labels) << std::endl;
  if (out_stream.bad()) {
    Failwith("Failure writing to " + out_path);
  }
  out_stream.close();
}

// ** TP Engine

void GPInstance::MakeTPEngine() {
  auto site_pattern = MakeSitePattern();
  std::string mmap_likelihood_path = mmap_file_path_.value() + ".tp_lik";
  std::string mmap_parsimony_path = mmap_file_path_.value() + ".tp_pars";
  tp_engine_ = std::make_unique<TPEngine>(GetDAG(), site_pattern, mmap_likelihood_path,
                                          mmap_parsimony_path);
}

TPEngine &GPInstance::GetTPEngine() {
  Assert(tp_engine_,
         "TPEngine not available. Call MakeTPEngine before accessing TPEngine.");
  return *tp_engine_;
}

void GPInstance::TPEngineSetChoiceMapByTakingFirst(const bool use_subsplit_method) {
  GetTPEngine().SetChoiceMapByTakingFirst(
      GetCurrentlyLoadedTrees(), GetDAG().BuildEdgeIndexer(), use_subsplit_method);
}

void GPInstance::TPEngineSetBranchLengthsByTakingFirst() {
  GetTPEngine().SetBranchLengthsByTakingFirst(GetCurrentlyLoadedTrees(),
                                              GetDAG().BuildEdgeIndexer());
}

// ** NNI Engine

void GPInstance::MakeNNIEngine() {
  nni_engine_ = std::make_unique<NNIEngine>(GetDAG(), nullptr, nullptr);
  if (gp_engine_ != nullptr) {
    GetNNIEngine().MakeGPEvalEngine(gp_engine_.get());
  }
  if (tp_engine_ != nullptr) {
    GetNNIEngine().MakeTPEvalEngine(tp_engine_.get());
  }
}

NNIEngine &GPInstance::GetNNIEngine() {
  Assert(nni_engine_,
         "NNIEngine not available. Call MakeNNIEngine before accessing NNIEngine.");
  return *nni_engine_;
}
