// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_instance.hpp"

#include "driver.hpp"
#include "gp_operation.hpp"
#include "numerical_utils.hpp"

#include <stdio.h>
#include <deque>

using namespace GPOperations;

void GPInstance::ReadFastaFile(std::string fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

void GPInstance::ReadNewickFile(std::string fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
}

void GPInstance::ReadNexusFile(std::string fname) {
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

void GPInstance::MakeEngine() {
  CheckSequencesAndTreesLoaded();
  ProcessLoadedTrees();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
  //std::cout << sbn_parameters_.size() << ", " << tree_collection_.TaxonCount() << std::endl;

  // SHJ: This is not quite correct. There can be more parameters involving
  // the trivial (fake) subsplits.
  engine_ = std::make_unique<GPEngine>(
      // To count GPCSSs, we add the usual suspects to the number of leaves (which are
      // the fake PCSS).
      site_pattern, sbn_parameters_.size() + tree_collection_.TaxonCount(),
      mmap_file_path_);
}

void GPInstance::MakeGPEngine() {
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());

  // Get number of parameters involving fake subsplits.
  size_t n_fake_subsplit_params = 0;
  for (size_t i = 0; i < alignment_.SequenceCount(); i++) {
    n_fake_subsplit_params += dag_nodes_[i]->GetRootwardRotated().size();
    n_fake_subsplit_params += dag_nodes_[i]->GetRootwardSorted().size();
  }

  std::cout << "Num fake split params: " << n_fake_subsplit_params << std::endl;

  //size_t n_leaves = alignment_.SequenceCount();
  size_t n_nodes = dag_nodes_.size();
  // The number of PLV needed is as follows:
  // 2 for the fake subsplit s: one for the observation phat(s) and one for p(s).
  // 2 for the pseudo root: one to combine the results from rootward pass, and
  // one for the stationary.
  // 6 for the other nodes: phat(s), phat(s.rotated), p(s) and r(s),
  // r(s.rotated), and rhat(s).
  //size_t n_plvs = 2 * n_leaves + 6 * (n_nodes - 1 - n_leaves) + 2;
  size_t n_plvs = 6 * n_nodes;
  engine_ = std::make_unique<GPEngine>(site_pattern,
                                       n_plvs,
                                       sbn_parameters_.size() + n_fake_subsplit_params,
                                       mmap_file_path_);
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

void GPInstance::ClearTreeCollectionAssociatedState() {
  sbn_parameters_.resize(0);
  rootsplits_.clear();
  indexer_.clear();
  index_to_child_.clear();
  parent_to_range_.clear();
}

void GPInstance::ProcessLoadedTrees() {
  size_t index = 0;
  ClearTreeCollectionAssociatedState();
  auto topology_counter = tree_collection_.TopologyCounter();
  // Start by adding the rootsplits.
  for (const auto &iter : RootedSBNMaps::RootsplitCounterOf(topology_counter)) {
    SafeInsert(indexer_, iter.first, index);
    rootsplits_.push_back(iter.first);
    index++;
  }
  // Now add the PCSSs.
  for (const auto &[parent, child_counter] :
       RootedSBNMaps::PCSSCounterOf(topology_counter)) {
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
}

void GPInstance::CreateAndInsertNode(const Bitset &subsplit)
{
  if (!subsplit_to_index_.count(subsplit)) {
    size_t id = dag_nodes_.size();
    subsplit_to_index_[subsplit] = id;
    dag_nodes_.push_back(std::make_shared<DAGNode>(id, subsplit));
    std::cout << id << ":" << subsplit.SubsplitToString() << std::endl;
  }
}

void GPInstance::AddChildrenSubsplits(const Bitset &subsplit,
                                      std::deque<Bitset> &q)
{
  // Retrieve children subsplits, add it to the queue to be processed.
  auto children_subsplits = GetChildrenSubsplits(subsplit);
  for (auto child_subsplit : children_subsplits) {
    q.push_back(child_subsplit);
  }
}

void GPInstance::ConnectNodes(size_t idx,
                              bool rotated)
{
  auto node = dag_nodes_[idx];
  // Retrieve children subsplits, set edge relation.
  Bitset subsplit = rotated ? node->GetBitset().RotateSubsplit() :
                              node->GetBitset();
  EdgeType leaf_edge_type = rotated ? EdgeType::LEAFWARD_ROTATED :
                                      EdgeType::LEAFWARD_SORTED;
  EdgeType root_edge_type = rotated ? EdgeType::ROOTWARD_ROTATED :
                                      EdgeType::ROOTWARD_SORTED;

  auto children = GetChildrenSubsplits(subsplit, true);
  for (auto child_subsplit : children) {
    auto child_node = dag_nodes_[subsplit_to_index_[child_subsplit]];
    //std::cout << child_node->Id() << std::endl;
    node->AddNeighbor(leaf_edge_type, child_node->Id());
    child_node->AddNeighbor(root_edge_type, node->Id());
  }

  //std::cout << node->ToString() << std::endl;
}

std::vector<Bitset> GPInstance::GetChildrenSubsplits(const Bitset &subsplit,
                                                     bool include_fake_subsplits)
{
  std::vector<Bitset> children_subsplits;
  
//  if (!subsplit.SplitChunk(1).Any()) {
//    return children_subsplits;
//  }

  if (parent_to_range_.count(subsplit)) {
    auto range = parent_to_range_.at(subsplit);
    for (auto idx = range.first; idx < range.second; idx++) {
      auto child_subsplit = index_to_child_.at(idx);
      children_subsplits.push_back(child_subsplit);
    }
  }
  else {
    // In the case where second chunk of the subsplit is a trivial subsplit,
    // it will not map to any value (parent_to_range_[subsplit] doesn't exist).
    // But we still need to create and connect to fake subsplits in the DAG.
    // TODO: Assertion that only one bit is 1 for subsplit.SplitChunk(1)
    if (include_fake_subsplits && subsplit.SplitChunk(0).Any()) {
      // The trivial subsplit is of the second chunk of subsplit.
      // Prepend it by 0's.
      Bitset zero(subsplit.size()/2);
      Bitset fake_subsplit = zero + subsplit.SplitChunk(1);
      //std::cout << trivial_subsplit.SubsplitToString() << std::endl;
      children_subsplits.push_back(fake_subsplit);
    }
  }
  return children_subsplits;
}

void GPInstance::BuildNodes()
{
  // We will create fake subsplits and insert to dag_nodes_.
  size_t n_taxa = this->alignment_.SequenceCount();
  Bitset zero(n_taxa);
  for (size_t i = 0; i < n_taxa; i++) {
    Bitset fake(n_taxa);
    fake.set(i);
    Bitset fake_subsplit = zero + fake;
    CreateAndInsertNode(fake_subsplit);
  }
  
  std::deque<Bitset> q;

  // We fill the next entries of dag_nodes_ using the root subsplits. subsplits.
  // And, populate the queue with children of rootsplits.
  for (auto rootsplit : rootsplits_) {
    auto subsplit = rootsplit + ~rootsplit;
    CreateAndInsertNode(subsplit);
    AddChildrenSubsplits(subsplit, q);
    AddChildrenSubsplits(subsplit.RotateSubsplit(), q);
  }
  
  // Fill the rest of dag_nodes_ with other subsplits.
  std::unordered_set<Bitset> visisted_subsplits;
  while (!q.empty()) {
    auto subsplit = q.front();
    q.pop_front();
    CreateAndInsertNode(subsplit);
    AddChildrenSubsplits(subsplit, q);
    AddChildrenSubsplits(subsplit.RotateSubsplit(), q);
  }
  
  // We will create a pseudoroot so that we can pass stationary distribution.
  CreateAndInsertNode(Bitset(2*rootsplits_.at(0).size(), true));
}

void GPInstance::BuildEdges()
{
  auto pseudoroot = dag_nodes_[dag_nodes_.size()-1];
  size_t n_fake_subsplits = this->alignment_.SequenceCount();
  for (size_t i = n_fake_subsplits; i < dag_nodes_.size() - 1; i++) {
    ConnectNodes(i, false);
    ConnectNodes(i, true);
    
    if (dag_nodes_[i]->IsRoot()) {
      pseudoroot->AddNeighbor(EdgeType::LEAFWARD_SORTED, dag_nodes_[i]->Id());
      dag_nodes_[i]->AddNeighbor(EdgeType::ROOTWARD_SORTED, pseudoroot->Id());
    }
  }
}

bool IsFakeSubsplit(const Bitset &subsplit) {
  return !subsplit.SplitChunk(0).Any();
}

void GPInstance::ConstructDAG()
{
  BuildNodes();
  BuildEdges();
}

void GPInstance::PrintDAG() {
  for (size_t i = 0; i < dag_nodes_.size(); i++) {
    std::cout << dag_nodes_[i]->ToString() << std::endl;
  }
}

void RootwardDepthFirst(size_t id,
                        std::vector<std::shared_ptr<DAGNode>> &dag_nodes,
                        std::vector<size_t> &visit_order,
                        std::unordered_set<size_t> &visited_nodes)
{
  visited_nodes.insert(id);
  for (size_t child_id : dag_nodes.at(id)->GetRootwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : dag_nodes.at(id)->GetRootwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

void LeafwardDepthFirst(size_t id,
                        std::vector<std::shared_ptr<DAGNode>> &dag_nodes,
                        std::vector<size_t> &visit_order,
                        std::unordered_set<size_t> &visited_nodes)
{
  visited_nodes.insert(id);
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

std::vector<size_t> GPInstance::LeafwardPassTraversal()
{
  size_t n_taxa = alignment_.SequenceCount();
//  Assert(leaf_idx < n_taxa, std::to_string(leaf_idx) + " >= num taxa.");
  std::vector<size_t> visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (size_t leaf_idx = 0; leaf_idx < n_taxa; leaf_idx++) {
    RootwardDepthFirst(leaf_idx,
                       dag_nodes_,
                       visit_order,
                       visited_nodes);
  }
  std::cout << "Leafward traversal node order:\n";
  for (size_t i = 0; i < visit_order.size(); i++) {
    std::cout << visit_order[i] << std::endl;
  }
  std::cout << "\n";
  return visit_order;
}

std::vector<size_t> GPInstance::RootwardPassTraversal()
{
  std::vector<size_t> visit_order;
  std::unordered_set<size_t> visited_nodes;
  LeafwardDepthFirst(dag_nodes_.size() - 1,
                     dag_nodes_,
                     visit_order,
                     visited_nodes);
  std::cout << "Rootward traversal node order:\n";
  for (size_t i = 0; i < visit_order.size(); i++) {
    std::cout << visit_order[i] << std::endl;
  }
  std::cout << "\n";
  return visit_order;
}

void GPInstance::BuildGPEngineIndexer(BitsetSizeMap &indexer,
                                      BitsetSizePairMap &subsplit2range)
{
  auto pseudoroot = dag_nodes_.back();
  size_t n_taxa = alignment_.SequenceCount();
  //SafeInsert(subsplit2range, pseudoroot->GetBitset(), {0, rootsplits_.size()});
  //size_t idx = rootsplits_.size();
  size_t idx = 0;
  for (size_t i = n_taxa; i < dag_nodes_.size(); i++) {
    auto node = dag_nodes_[i];
    std::cout << node->ToString() << std::endl;
    std::cout << node->GetBitset().SubsplitToString() << std::endl;
    auto n_child = node->GetLeafwardSorted().size();
    if (n_child > 0) {
      SafeInsert(subsplit2range, node->GetBitset(), {idx, idx + n_child});
      for (size_t j = 0; j < n_child; j++) {
        auto child = dag_nodes_[node->GetLeafwardSorted()[j]];
        indexer[node->GetBitset() + child->GetBitset()] = idx;
        idx++;
      }
    }
    n_child = node->GetLeafwardRotated().size();
    if (n_child > 0) {
      SafeInsert(subsplit2range, node->GetBitset().RotateSubsplit(), {idx, idx + n_child});
      for (size_t j = 0; j < n_child; j++) {
        auto child = dag_nodes_[node->GetLeafwardRotated()[j]];
        indexer[node->GetBitset().RotateSubsplit() + child->GetBitset()] = idx;
        //indexer[node->GetBitset() + child->GetBitset()] = idx;
        idx++;
      }
    }
  }
}

void GPInstance::RootwardPass(std::vector<size_t> visit_order,
                              bool optimize)
{
  // Perform first rootward pass. No optimization.
  size_t n_nodes = dag_nodes_.size();
  size_t pcsp_idx = 0;
  GPOperationVector operations;
  for (size_t src_idx : visit_order) {
    auto node = dag_nodes_[src_idx];
    if (node->IsLeaf()) {
    } else if (node->IsRoot()) {
//      operations.push_back(CopyPLV{src_idx, n_nodes + src_idx});
//      operations.push_back(CopyPLV{2*n_nodes + src_idx, n_nodes + src_idx});
    } else {
      // Combine phat(s), phat(s_tilde) and place it into p(s).
      operations.push_back(Multiply{src_idx, n_nodes + src_idx, 2*n_nodes + src_idx});

      // Update the likelihood.
      for (size_t child_idx : node->GetLeafwardSorted()) {
        auto child_node = dag_nodes_[child_idx];
        pcsp_idx = gp_engine_indexer_[node->GetBitset() + child_node->GetBitset()];
        // Update likelihood l(s|t).
        operations.push_back(Likelihood{pcsp_idx, child_idx, 4*n_nodes + src_idx});
      }
      for (size_t child_idx : node->GetLeafwardRotated()) {
        auto child_node = dag_nodes_[child_idx];
        pcsp_idx = gp_engine_indexer_[node->GetBitset().RotateSubsplit() + child_node->GetBitset()];
        // Update likelihood l(s|t).
        operations.push_back(Likelihood{pcsp_idx, child_idx, 5*n_nodes + src_idx});
      }
    }

    for (size_t dest_idx : node->GetRootwardSorted()) {
      auto parent_node = dag_nodes_[dest_idx];
      pcsp_idx = gp_engine_indexer_[parent_node->GetBitset() + node->GetBitset()];
      if (!parent_node->IsRoot()) {
        operations.push_back(WeightedSumAccumulate{n_nodes + dest_idx, pcsp_idx,
          pcsp_idx, src_idx});
      } else {
        operations.push_back(WeightedSumAccumulateStationary{pcsp_idx,
          4*n_nodes + src_idx, src_idx});
      }
    }
    for (size_t dest_idx : node->GetRootwardRotated()) {
      auto parent_node = dag_nodes_[dest_idx];
      pcsp_idx = gp_engine_indexer_[parent_node->GetBitset().RotateSubsplit() + node->GetBitset()];
      operations.push_back(WeightedSumAccumulate{2*n_nodes + dest_idx, pcsp_idx,
        pcsp_idx, src_idx});
    }
    
//    if (optimize) {
//      if (subsplit2range_.count(node->GetBitset())) {
//        auto param_range = subsplit2range_[node->GetBitset()];
//        operations.push_back(UpdateSBNProbabilities{param_range.first, param_range.second});
//      }
//      if (subsplit2range_.count(node->GetBitset().RotateSubsplit())) {
//        auto param_range = subsplit2range_[node->GetBitset().RotateSubsplit()];
//        operations.push_back(UpdateSBNProbabilities{param_range.first, param_range.second});
//      }
//    }
  }

  GetEngine()->ProcessOperations(operations);
}

void GPInstance::LeafwardPass(std::vector<size_t> visit_order,
                              bool optimize)
{
  size_t n_nodes = dag_nodes_.size();
  size_t pcsp_idx;
  GPOperationVector operations;
  for (size_t src_idx : visit_order) {
    auto node = dag_nodes_[src_idx];
    if (!node->IsRoot()) {
      // Fill r(s) = phat(s_tilde) \circ rhat(s).
      // Recall:
      // phat(s_tilde) is stored in [2*n_nodes, 3*n_nodes),
      // rhat(s) is stored in [3*n_nodes, 4*n_nodes), and
      // r(s) is stored in [4*n_nodes, 5*n_nodes).
      operations.push_back(Multiply{4*n_nodes + src_idx, 3*n_nodes + src_idx, 2*n_nodes + src_idx});
      // Fill r(s_tilde) = phat(s) \circ rhat(s).
      // Recall:
      // phat(s) is stored in [n_nodes, 2*n_nodes),
      // rhat(s) is stored in [3*n_nodes, 4*n_nodes), and
      // r(s_tilde) is stored in [5*n_nodes, 6*n_nodes).
      operations.push_back(Multiply{5*n_nodes + src_idx, 3*n_nodes + src_idx, n_nodes + src_idx});

      for (size_t dest_idx : node->GetLeafwardSorted()) {
        auto child_node = dag_nodes_[dest_idx];
        pcsp_idx = gp_engine_indexer_[node->GetBitset() + child_node->GetBitset()];
        // Update likelihood l(s|t).
        operations.push_back(Likelihood{pcsp_idx, dest_idx, 4*n_nodes + src_idx});
        // Update rhat(s) for each child s.
        operations.push_back(WeightedSumAccumulate{3*n_nodes + dest_idx, pcsp_idx,
          pcsp_idx, 4*n_nodes + src_idx});
      }
      for (size_t dest_idx : node->GetLeafwardRotated()) {
        auto child_node = dag_nodes_[dest_idx];
        pcsp_idx = gp_engine_indexer_[node->GetBitset().RotateSubsplit() + child_node->GetBitset()];
        // Update likelihood l(s|t).
        operations.push_back(Likelihood{pcsp_idx, dest_idx, 5*n_nodes + src_idx});
        // Update rhat(s) for each child s_tilde.
        operations.push_back(WeightedSumAccumulate{3*n_nodes + dest_idx, pcsp_idx,
          pcsp_idx, 5*n_nodes + src_idx});
      }
    } else {
      // For pseudoroot, we just need to compute the likelihood.
      // q(s) \pi p(dest_idx)
      for (size_t dest_idx : node->GetLeafwardSorted()) {
        auto child_node = dag_nodes_[dest_idx];
        pcsp_idx = gp_engine_indexer_[node->GetBitset() + child_node->GetBitset()];
        operations.push_back(Likelihood{pcsp_idx, dest_idx, 4*n_nodes + dest_idx});
      }
    }

    if (optimize) {
      if (subsplit2range_.count(node->GetBitset())) {
        auto param_range = subsplit2range_[node->GetBitset()];
        operations.push_back(UpdateSBNProbabilities{param_range.first, param_range.second});
      }
      if (!node->IsRoot()) {
        if (subsplit2range_.count(node->GetBitset().RotateSubsplit())) {
          auto param_range = subsplit2range_[node->GetBitset().RotateSubsplit()];
          operations.push_back(UpdateSBNProbabilities{param_range.first, param_range.second});
        }
      }
    }
  }

  GetEngine()->ProcessOperations(operations);
}

void GPInstance::InitializeLikelihoodEM()
{
  size_t n_params = engine_->GetBranchLengths().size();

  // Initialize branch lengths to 1.
  EigenVectorXd branch_lengths(n_params);
  branch_lengths = EigenVectorXd::Ones(n_params);
  GetEngine()->SetBranchLengths(branch_lengths);

  // TODO: Initialize SBN parameters using simple average.
  // For now, use uniform.
  EigenVectorXd q = EigenVectorXd::Ones(n_params);
  for (auto it = subsplit2range_.begin(); it != subsplit2range_.end(); ++it) {
    auto range = subsplit2range_[it->first];
    auto num_child_subsplits = range.second - range.first;
    double val = 1./num_child_subsplits;
    q.segment(range.first, num_child_subsplits).array() = val;
  }
  GetEngine()->SetSBNParameters(q);
  std::cout << q << std::endl;

  // Set stationary for the pseudoroot.
  std::vector<GPOperation> operations;
  size_t n_nodes = dag_nodes_.size();
  
  // For each leaf node, we have p(s) populated by GPEngine.
  // For completeness, we fill phat(s_tilde) = phat(s) and p(s) = phat(s).
  size_t n_taxa = this->alignment_.SequenceCount();
  for (size_t src_idx = 0; src_idx < n_taxa; src_idx++) {
    auto node = dag_nodes_[src_idx];
    operations.push_back(CopyPLV{n_nodes + src_idx, src_idx});
    operations.push_back(CopyPLV{2*n_nodes + src_idx, src_idx});
  }
  // For the rootsplits s, we want to fill rhat(s) = stationary.
  for (auto rootsplit_idx : dag_nodes_[n_nodes - 1]->GetLeafwardSorted()) {
    auto node = dag_nodes_[rootsplit_idx];
    operations.push_back(SetToStationaryDistribution{3*n_nodes + rootsplit_idx});
    operations.push_back(SetToStationaryDistribution{4*n_nodes + rootsplit_idx});
    operations.push_back(SetToStationaryDistribution{5*n_nodes + rootsplit_idx});
  }
    
  GetEngine()->ProcessOperations(operations);
}

void GPInstance::SetZero()
{
  GPOperationVector operations;
  auto n_nodes = dag_nodes_.size();
  auto n_taxa = alignment_.SequenceCount();
  for (size_t i = 0; i < n_taxa; i++) {
    for (size_t j = 3; j < 6; j++) {
      operations.push_back(Zero{j * n_nodes + i});
    }
  }
  for (size_t i = n_taxa; i < n_nodes-1; i++) {
    for (size_t j = 0; j < 6; j++) {
      operations.push_back(Zero{j * n_nodes + i});
    }
  }
  for (auto rootsplit_idx : dag_nodes_[n_nodes - 1]->GetLeafwardSorted()) {
    auto node = dag_nodes_[rootsplit_idx];
    operations.push_back(SetToStationaryDistribution{3*n_nodes + rootsplit_idx});
    operations.push_back(SetToStationaryDistribution{4*n_nodes + rootsplit_idx});
    operations.push_back(SetToStationaryDistribution{5*n_nodes + rootsplit_idx});
  }
  GetEngine()->ProcessOperations(operations);
}


void GPInstance::TrainLikelihoodEM(double tol, size_t max_iter) {
  CheckSequencesAndTreesLoaded();
  ProcessLoadedTrees();

  ConstructDAG();
  //PrintDAG();

  // Construct GPEngine.
  MakeGPEngine();

  // Build indexer for q_, branch_lengths_, and log_likelihoods_ in GPEngine.
  BuildGPEngineIndexer(gp_engine_indexer_, subsplit2range_);

  InitializeLikelihoodEM();

  auto rootward_order = RootwardPassTraversal();
  auto leafward_order = LeafwardPassTraversal();

  // Perform 1st rootward pass to establish likelihood for the current parameters.
  RootwardPass(rootward_order, false);
  // Reset marginal likelihood.
  GetEngine()->ResetLogMarginalLikelihood();
  double curr_log_lik = DOUBLE_NEG_INF;

  // EM Loop: LeafwardPass, RootwardPass.
  for (size_t i = 0; i < max_iter; i++) {
    std::cout << "EM iter: " << i << std::endl;
    LeafwardPass(leafward_order, true);
    for (size_t i = 3; i < 6; i++) {
      for (size_t j = 0; j < dag_nodes_.size(); j++) {
        GetEngine()->PrintPLV(i * dag_nodes_.size() + j);
      }
    }
    SetZero();
    GetEngine()->ResetLogMarginalLikelihood();
    RootwardPass(rootward_order, true);
    double log_lik = GetEngine()->GetLogMarginalLikelihood();
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < dag_nodes_.size(); j++) {
        GetEngine()->PrintPLV(i * dag_nodes_.size() + j);
      }
    }
    std::cout << "LogLik: " << log_lik << std::endl;

    // Get loglikelihood and check for convergence.
    if (abs(log_lik - curr_log_lik) < tol) {
      break;
    }
    curr_log_lik = log_lik;
  }
}
