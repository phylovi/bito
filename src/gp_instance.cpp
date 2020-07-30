// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_instance.hpp"

#include <chrono>
#include <iomanip>
#include <string>

#include "driver.hpp"
#include "gp_operation.hpp"
#include "numerical_utils.hpp"

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
  std::cout << dag_.NodeCount() << " DAG nodes representing " << dag_.TreeCount()
            << " trees.\n";
  std::cout << dag_.GeneralizedPCSPCount() << " continuous parameters.\n";
  if (engine_ == nullptr) {
    std::cout << "Engine has not been made.\n";
  } else {
    std::cout << "Engine available with " << GetEngine()->PLVByteCount() / 1e9
              << "G virtual memory.\n";
  }
  
  for (auto &trees : dag_.GenerateAllTrees()) {
    std::cout << trees << "\n";
  }
}

void GPInstance::ReadFastaFile(const std::string &fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

void GPInstance::ReadNewickFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
}

void GPInstance::ReadNexusFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
}

void GPInstance::CheckSequencesAndTreesLoaded() const {
  if (alignment_.SequenceCount() == 0) {
    Failwith(
        "Load an alignment into your GPInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
  if (tree_collection_.TreeCount() == 0) {
    Failwith(
        "Load some trees into your GPInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
}

void GPInstance::MakeEngine(double rescaling_threshold) {
  CheckSequencesAndTreesLoaded();
  ProcessLoadedTrees();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());

  dag_ = GPDAG(tree_collection_);
  engine_ = std::make_unique<GPEngine>(site_pattern, 6 * dag_.NodeCount(),
                                       dag_.GeneralizedPCSPCount(), mmap_file_path_,
                                       rescaling_threshold);
  InitializeGPEngine();
}

GPEngine *GPInstance::GetEngine() const {
  if (engine_ != nullptr) {
    return engine_.get();
  }
  // else
  Failwith(
      "Engine not available. Call MakeEngine to make an engine for phylogenetic "
      "likelihood computation.");
}

void GPInstance::ProcessOperations(const GPOperationVector &operations) {
  GetEngine()->ProcessOperations(operations);
}

void GPInstance::ClearTreeCollectionAssociatedState() {
  sbn_parameters_.resize(0);
  dag_ = GPDAG();
}

void GPInstance::ProcessLoadedTrees() {
  sbn_parameters_.resize(dag_.RootsplitAndPCSPCount());
  sbn_parameters_.setOnes();
}

void GPInstance::PrintDAG() { dag_.Print(); }
void GPInstance::PrintGPCSPIndexer() { dag_.PrintGPCSPIndexer(); }

void GPInstance::InitializeGPEngine() {
  GetEngine()->SetSBNParameters(dag_.BuildUniformPrior());
}

void GPInstance::PopulatePLVs() {
  ProcessOperations(dag_.SetRootwardZero());
  ProcessOperations(dag_.SetLeafwardZero());
  ProcessOperations(dag_.SetRhatToStationary());
  ProcessOperations(dag_.RootwardPass());
  ProcessOperations(dag_.LeafwardPass());
}

void GPInstance::ComputeLikelihoods() { ProcessOperations(dag_.ComputeLikelihoods()); }

void GPInstance::EstimateBranchLengths(double tol, size_t max_iter) {
  auto now = std::chrono::high_resolution_clock::now;
  auto t_start = now();
  std::cout << "Begin branch optimization\n";
  GPOperationVector branch_optimization_operations = dag_.BranchLengthOptimization();
  GPOperationVector marginal_lik_operations = dag_.MarginalLikelihood();

  GetEngine()->ResetLogMarginalLikelihood();
  std::cout << "Populating PLVs\n";
  PopulatePLVs();
  std::chrono::duration<double> warmup_duration = now() - t_start;
  t_start = now();
  std::cout << "Computing initial likelihood\n";
  ProcessOperations(marginal_lik_operations);
  double current_marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
  std::chrono::duration<double> initial_likelihood_duration = now() - t_start;
  t_start = now();

  for (size_t i = 0; i < max_iter; i++) {
    std::cout << "Iteration: " << (i + 1) << std::endl;
    ProcessOperations(branch_optimization_operations);
    GetEngine()->ResetLogMarginalLikelihood();
    ProcessOperations(marginal_lik_operations);
    double marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
    std::cout << "Current marginal log likelihood: ";
    std::cout << std::setprecision(9) << current_marginal_log_lik << std::endl;
    std::cout << "New marginal log likelihood: ";
    std::cout << std::setprecision(9) << marginal_log_lik << std::endl;
    if (marginal_log_lik < current_marginal_log_lik) {
      std::cout << "Marginal log likelihood decreased. Check branch optimization: "
                   "adjust step size.\n";
    }
    if (abs(current_marginal_log_lik - marginal_log_lik) < tol) {
      std::cout << "Converged.\n";
      break;
    }
    current_marginal_log_lik = marginal_log_lik;
  }
  std::chrono::duration<double> optimization_duration = now() - t_start;
  std::cout << "\n# Timing Report\n";
  std::cout << "warmup: " << warmup_duration.count() << "s\n";
  std::cout << "initial likelihood: " << initial_likelihood_duration.count() << "s\n";
  std::cout << "optimization: " << optimization_duration.count() << "s or "
            << optimization_duration.count() / 60 << "m\n";
}

void GPInstance::EstimateSBNParameters() {
  std::cout << "Begin SBN parameter optimization\n";
  GPOperationVector sbn_param_optimization_operations = dag_.OptimizeSBNParameters();
  GPOperationVector marginal_lik_operations = dag_.MarginalLikelihood();

  GetEngine()->ResetLogMarginalLikelihood();
  PopulatePLVs();
  ComputeLikelihoods();

  ProcessOperations(sbn_param_optimization_operations);
  ProcessOperations(marginal_lik_operations);
  double marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
  std::cout << std::setprecision(9) << marginal_log_lik << std::endl;
}
