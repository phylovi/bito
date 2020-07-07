// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_instance.hpp"

#include <iomanip>
#include <stdio.h>
#include <string>

#include "driver.hpp"
#include "gp_operation.hpp"
#include "numerical_utils.hpp"

using namespace GPOperations;

void GPInstance::PrintStatus() {
  const auto tree_count = tree_collection_.TreeCount();
  const auto taxon_count = tree_collection_.TaxonCount();
  if (tree_count) {
    std::cout << tree_count << " unique tree topologies loaded on " << taxon_count
              << " leaves.\n";
  } else {
    std::cout << "No trees loaded.\n";
  }
  std::cout << alignment_.Data().size() << " sequences loaded.\n";
}

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

  ConstructDAG();

  MakeGPEngine();
  BuildPCSPIndexer();

  rootward_order_ = RootwardPassTraversal();
  leafward_order_ = LeafwardPassTraversal();

  InitializeGPEngine();
}

void GPInstance::MakeGPEngine() {
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());

  // Get number of parameters involving fake subsplits.
  size_t n_fake_subsplit_params = 0;
  for (size_t i = 0; i < alignment_.SequenceCount(); i++) {
    n_fake_subsplit_params += dag_nodes_[i]->GetRootwardRotated().size();
    n_fake_subsplit_params += dag_nodes_[i]->GetRootwardSorted().size();
  }

  size_t node_count = dag_nodes_.size();
  size_t plv_count = 6 * node_count;
  engine_ = std::make_unique<GPEngine>(site_pattern,
                                       plv_count,
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
  }
}

void GPInstance::AddChildrenSubsplits(const Bitset &subsplit,
                                      std::deque<Bitset> &q,
                                      std::unordered_set<Bitset> &visited_subsplits)
{
  // Retrieve children subsplits, add it to the queue to be processed.
  auto children_subsplits = GetChildrenSubsplits(subsplit);
  for (auto child_subsplit : children_subsplits) {
    if (!visited_subsplits.count(child_subsplit)) {
      q.push_back(child_subsplit);
      visited_subsplits.insert(child_subsplit);
    }
  }
}

void GPInstance::ConnectNodes(size_t idx,
                              bool rotated)
{
  auto node = dag_nodes_[idx];
  // Retrieve children subsplits, set edge relation.
  Bitset subsplit = rotated ? node->GetBitset().RotateSubsplit() :
                              node->GetBitset();
  auto children = GetChildrenSubsplits(subsplit, true);
  for (auto child_subsplit : children) {
    auto child_node = dag_nodes_[subsplit_to_index_[child_subsplit]];
    if (rotated) {
      node->AddLeafwardRotated(child_node->Id());
      child_node->AddRootwardRotated(node->Id());
    } else {
      node->AddLeafwardSorted(child_node->Id());
      child_node->AddRootwardSorted(node->Id());
    }
  }
}

std::vector<Bitset> GPInstance::GetChildrenSubsplits(const Bitset &subsplit,
                                                     bool include_fake_subsplits)
{
  std::vector<Bitset> children_subsplits;

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
    // TODO: Assertion that exactly one bit is 1 for subsplit.SplitChunk(1)
    if (include_fake_subsplits && subsplit.SplitChunk(0).Any()) {
      // The fake subsplit corresponds to the second chunk of subsplit.
      // Prepend it by 0's.
      Bitset zero(subsplit.size()/2);
      Bitset fake_subsplit = zero + subsplit.SplitChunk(1);
      children_subsplits.push_back(fake_subsplit);
    }
  }
  return children_subsplits;
}

void GPInstance::BuildNodes()
{
  // We will create fake subsplits and insert to dag_nodes_.
  size_t taxon_count = this->alignment_.SequenceCount();
  Bitset zero(taxon_count);
  for (size_t i = 0; i < taxon_count; i++) {
    Bitset fake(taxon_count);
    fake.set(i);
    Bitset fake_subsplit = zero + fake;
    CreateAndInsertNode(fake_subsplit);
  }

  std::deque<Bitset> q;
  std::unordered_set<Bitset> visited_subsplits;

  // We fill the next entries of dag_nodes_ using the root subsplits. subsplits.
  // And, populate the queue with children of rootsplits.
  for (auto rootsplit : rootsplits_) {
    auto subsplit = rootsplit + ~rootsplit;
    CreateAndInsertNode(subsplit);
    AddChildrenSubsplits(subsplit, q, visited_subsplits);
    AddChildrenSubsplits(subsplit.RotateSubsplit(), q, visited_subsplits);
  }

  // Fill the rest of dag_nodes_ with other subsplits.
  while (!q.empty()) {
    auto subsplit = q.front();
    q.pop_front();
    CreateAndInsertNode(subsplit);
    AddChildrenSubsplits(subsplit, q, visited_subsplits);
    AddChildrenSubsplits(subsplit.RotateSubsplit(), q, visited_subsplits);
  }
}

void GPInstance::BuildEdges()
{
  size_t taxon_count = this->alignment_.SequenceCount();
  for (size_t i = taxon_count; i < dag_nodes_.size(); i++) {
    ConnectNodes(i, false);
    ConnectNodes(i, true);
  }
}

void GPInstance::ConstructDAG()
{
  BuildNodes();
  BuildEdges();
//  PrintDAG();
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
  size_t taxon_count = alignment_.SequenceCount();
//  Assert(leaf_idx < taxon_count, std::to_string(leaf_idx) + " >= num taxa.");
  std::vector<size_t> visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (size_t leaf_idx = 0; leaf_idx < taxon_count; leaf_idx++) {
    RootwardDepthFirst(leaf_idx,
                       dag_nodes_,
                       visit_order,
                       visited_nodes);
  }
  
  return visit_order;
}

std::vector<size_t> GPInstance::RootwardPassTraversal()
{
  std::vector<size_t> visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (auto rootsplit : rootsplits_) {
    size_t root_idx = subsplit_to_index_[rootsplit + ~rootsplit];
    LeafwardDepthFirst(root_idx,
                       dag_nodes_,
                       visit_order,
                       visited_nodes);
  }
  return visit_order;
}

void GPInstance::BuildPCSPIndexer()
{
  size_t taxon_count = alignment_.SequenceCount();
  size_t idx = 0;
  for (auto rootsplit : rootsplits_) {
    SafeInsert(pcsp_indexer_, rootsplit + ~rootsplit, idx);
    idx++;
  }

  for (size_t i = taxon_count; i < dag_nodes_.size(); i++) {
    auto node = dag_nodes_[i];
    auto n_child = node->GetLeafwardSorted().size();
    if (n_child > 0) {
      SafeInsert(subsplit2range_, node->GetBitset(), {idx, idx + n_child});
      for (size_t j = 0; j < n_child; j++) {
        auto child = dag_nodes_[node->GetLeafwardSorted()[j]];
        SafeInsert(pcsp_indexer_, node->GetBitset() + child->GetBitset(), idx);
        idx++;
      }
    }
    n_child = node->GetLeafwardRotated().size();
    if (n_child > 0) {
      SafeInsert(subsplit2range_,
                 node->GetBitset().RotateSubsplit(),
                 {idx, idx + n_child});
      for (size_t j = 0; j < n_child; j++) {
        auto child = dag_nodes_[node->GetLeafwardRotated()[j]];
        SafeInsert(pcsp_indexer_,
                   node->GetBitset().RotateSubsplit() + child->GetBitset(),
                   idx);
        idx++;
      }
    }
  }
}

void GPInstance::InitializeGPEngine()
{
  size_t n_params = engine_->GetBranchLengths().size();

  // Initialize branch lengths to 1.
  EigenVectorXd branch_lengths(n_params);
  branch_lengths = EigenVectorXd::Ones(n_params);
  GetEngine()->SetBranchLengths(branch_lengths);

  // For now, use uniform.
  EigenVectorXd q = EigenVectorXd::Ones(n_params);
  q.segment(0, rootsplits_.size()).array() = 1./rootsplits_.size();
  for (auto it = subsplit2range_.begin(); it != subsplit2range_.end(); ++it) {
    auto range = subsplit2range_[it->first];
    auto num_child_subsplits = range.second - range.first;
    double val = 1./num_child_subsplits;
    q.segment(range.first, num_child_subsplits).array() = val;
  }
  GetEngine()->SetSBNParameters(q);

  // Set rhat(s) = stationary for the rootsplits s.
  std::vector<GPOperation> operations;
  for (size_t i = 0; i < rootsplits_.size(); i++) {
    auto rootsplit = rootsplits_[i];
    auto root_idx = subsplit_to_index_[rootsplit + ~rootsplit];
    operations.push_back(SetToStationaryDistribution{
      GetPlvIndex(PlvType::R_HAT, dag_nodes_.size(), root_idx),
      i
    });
  }

  GetEngine()->ProcessOperations(operations);
}

size_t GetPlvIndex(PlvType plv_type,
                   size_t node_count,
                   size_t src_idx)
{
  switch (plv_type) {
    case PlvType::P:
      return src_idx;
    case PlvType::P_HAT:
      return node_count + src_idx;
    case PlvType::P_HAT_TILDE:
      return 2*node_count + src_idx;
    case PlvType::R_HAT:
      return 3*node_count + src_idx;
    case PlvType::R:
      return 4*node_count + src_idx;
    case PlvType::R_TILDE:
      return 5*node_count + src_idx;
    default:
      Failwith("Invalid PLV index requested.");
  }
}

void GPInstance::AddRootwardWeightedSumAccumulateOperations(
    std::shared_ptr<DAGNode> node,
    bool rotated,
    GPOperationVector &operations)
{
  std::vector<size_t> child_idxs = rotated ? node->GetLeafwardRotated() :
                                             node->GetLeafwardSorted();
  PlvType plv_type = rotated ? PlvType::P_HAT_TILDE : PlvType::P_HAT;
  auto parent_subsplit = rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  for (size_t child_idx : child_idxs) {
    auto child_node = dag_nodes_[child_idx];
    auto child_subsplit = child_node->GetBitset();
    auto pcsp = parent_subsplit + child_subsplit;
    if (!pcsp_indexer_.count(pcsp)) {
      Failwith("Non-existent PCSP index.");
    }
    auto pcsp_idx = pcsp_indexer_[pcsp];

    operations.push_back(WeightedSumAccumulate{
      GetPlvIndex(plv_type, dag_nodes_.size(), node->Id()),
      pcsp_idx,
      GetPlvIndex(PlvType::P, dag_nodes_.size(), child_idx)});
  }
}

void GPInstance::AddLeafwardWeightedSumAccumulateOperations(
     std::shared_ptr<DAGNode> node,
     GPOperationVector &operations)
{
  auto subsplit = node->GetBitset();
  for (size_t parent_idx : node->GetRootwardSorted()) {
    auto parent_node = dag_nodes_[parent_idx];
    auto parent_subsplit = parent_node->GetBitset();
    auto pcsp_idx = pcsp_indexer_[parent_subsplit + subsplit];

    operations.push_back(WeightedSumAccumulate{
      GetPlvIndex(PlvType::R_HAT, dag_nodes_.size(), node->Id()),
      pcsp_idx,
      GetPlvIndex(PlvType::R, dag_nodes_.size(), parent_node->Id())});
  }
  for (size_t parent_idx : node->GetRootwardRotated()) {
    auto parent_node = dag_nodes_[parent_idx];
    auto parent_subsplit = parent_node->GetBitset().RotateSubsplit();
    auto pcsp_idx = pcsp_indexer_[parent_subsplit + subsplit];

    operations.push_back(WeightedSumAccumulate{
      GetPlvIndex(PlvType::R_HAT, dag_nodes_.size(), node->Id()),
      pcsp_idx,
      GetPlvIndex(PlvType::R_TILDE, dag_nodes_.size(), parent_node->Id())});
  }
}

void GPInstance::OptimizeSBNParameters(const Bitset &subsplit,
                                       GPOperationVector &operations)
{
  if (subsplit2range_.count(subsplit)) {
    auto param_range = subsplit2range_[subsplit];
    if (param_range.second - param_range.first > 1) {
      operations.push_back(UpdateSBNProbabilities{param_range.first, param_range.second});
    }
  }
}

size_t GPInstance::GetPCSPIndex(size_t parent_node_idx,
                                size_t child_node_idx,
                                bool rotated) {
  auto parent_node = dag_nodes_[parent_node_idx];
  auto child_node = dag_nodes_[child_node_idx];
  auto pcsp = rotated ?
        parent_node->GetBitset().RotateSubsplit() + child_node->GetBitset() :
        parent_node->GetBitset() + child_node->GetBitset();
  return pcsp_indexer_[pcsp];
}

void GPInstance::RootwardPass(std::vector<size_t> visit_order)
{
  // Perform first rootward pass. No optimization.
  GPOperationVector operations;
  for (size_t node_idx : visit_order) {
    auto node = dag_nodes_[node_idx];
    if (node->IsLeaf()) {
      continue;
    }

    // Perform the following operations:
    // 1. WeightedSumAccummulate to get phat(s), phat(s_tilde).
    // 2. Multiply to get p(s) = phat(s) \circ phat(s_tilde).
    AddRootwardWeightedSumAccumulateOperations(node, false, operations);
    AddRootwardWeightedSumAccumulateOperations(node, true, operations);
    operations.push_back(Multiply{node_idx,
      GetPlvIndex(PlvType::P_HAT, dag_nodes_.size(), node_idx),
      GetPlvIndex(PlvType::P_HAT_TILDE, dag_nodes_.size(), node_idx)
    });
  }

  GetEngine()->ProcessOperations(operations);
}

void GPInstance::LeafwardPass(std::vector<size_t> visit_order)
{
  GPOperationVector operations;
  for (size_t node_idx : visit_order) {
    auto node = dag_nodes_[node_idx];

    // Perform the following operations:
    // 1. WeightedSumAccumulate: rhat(s) += \sum_t q(s|t) P'(s|t) r(t)
    // 2. Multiply: r(s) = rhat(s) \circ phat(s_tilde).
    // 3. Multiply: r(s_tilde) = rhat(s) \circ phat(s).
    AddLeafwardWeightedSumAccumulateOperations(node, operations);
    operations.push_back(Multiply{
      GetPlvIndex(PlvType::R, dag_nodes_.size(), node_idx),
      GetPlvIndex(PlvType::R_HAT, dag_nodes_.size(), node_idx),
      GetPlvIndex(PlvType::P_HAT_TILDE, dag_nodes_.size(), node_idx)});
    operations.push_back(Multiply{
      GetPlvIndex(PlvType::R_TILDE, dag_nodes_.size(), node_idx),
      GetPlvIndex(PlvType::R_HAT, dag_nodes_.size(), node_idx),
      GetPlvIndex(PlvType::P_HAT, dag_nodes_.size(), node_idx)});
  }
  GetEngine()->ProcessOperations(operations);
}

void GPInstance::SetRootwardZero()
{
  GPOperationVector operations;
  auto node_count = dag_nodes_.size();
  auto taxon_count = alignment_.SequenceCount();
  for (size_t i = taxon_count; i < node_count; i++) {
    operations.push_back(Zero{GetPlvIndex(PlvType::P, dag_nodes_.size(), i)});
    operations.push_back(Zero{GetPlvIndex(PlvType::P_HAT, dag_nodes_.size(), i)});
    operations.push_back(Zero{GetPlvIndex(PlvType::P_HAT_TILDE, dag_nodes_.size(), i)});
  }
  GetEngine()->ProcessOperations(operations);
}

void GPInstance::SetLeafwardZero()
{
  GPOperationVector operations;
  auto node_count = dag_nodes_.size();
  for (size_t i = 0; i < node_count; i++) {
    operations.push_back(Zero{GetPlvIndex(PlvType::R_HAT, dag_nodes_.size(), i)});
    operations.push_back(Zero{GetPlvIndex(PlvType::R, dag_nodes_.size(), i)});
    operations.push_back(Zero{GetPlvIndex(PlvType::R_TILDE, dag_nodes_.size(), i)});
  }
  for (auto rootsplit : rootsplits_) {
    size_t root_idx = subsplit_to_index_[rootsplit + ~rootsplit];
    operations.push_back(SetToStationaryDistribution{GetPlvIndex(PlvType::R_HAT,
                                                                 dag_nodes_.size(),
                                                                 root_idx)});
  }
  GetEngine()->ProcessOperations(operations);
}

void GPInstance::ScheduleBranchLengthOptimization(size_t node_id,
                                      std::unordered_set<size_t> &visited_nodes,
                                      GPOperationVector &operations) {
  visited_nodes.insert(node_id);
  auto node = dag_nodes_[node_id];

  if (!node->IsRoot()) {
    // Compute R_HAT(s) = \sum_{t : s < t} q(s|t) P(s|t) r(t).
    operations.push_back(Zero{GetPlvIndex(R_HAT, dag_nodes_.size(), node_id)});
    for (size_t parent_id : node->GetRootwardSorted()) {
      auto parent_node = dag_nodes_[parent_id];
      size_t pcsp_idx = pcsp_indexer_[parent_node->GetBitset() + node->GetBitset()];
      operations.push_back(WeightedSumAccumulate{
        GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
        pcsp_idx,
        GetPlvIndex(R, dag_nodes_.size(), parent_id)
      });
    }
    for (size_t parent_id : node->GetRootwardRotated()) {
      auto parent_node = dag_nodes_[parent_id];
      size_t pcsp_idx = pcsp_indexer_[parent_node->GetBitset().RotateSubsplit() + node->GetBitset()];
      operations.push_back(WeightedSumAccumulate{
        GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
        pcsp_idx,
        GetPlvIndex(R_TILDE, dag_nodes_.size(), parent_id)
      });
    }

    // Update r(s) and r_tilde(s)
    operations.push_back(Multiply{
      GetPlvIndex(R, dag_nodes_.size(), node_id),
      GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
      GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id)
    });
    operations.push_back(Multiply{
      GetPlvIndex(R_TILDE, dag_nodes_.size(), node_id),
      GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
      GetPlvIndex(P_HAT, dag_nodes_.size(), node_id)
    });
  }

  if (node->IsLeaf())
    return;

  operations.push_back(Zero{GetPlvIndex(P_HAT, dag_nodes_.size(), node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleBranchLengthOptimization(child_id, visited_nodes, operations);
    }

    auto child_node = dag_nodes_[child_id];
    size_t pcsp_idx = pcsp_indexer_[node->GetBitset() + child_node->GetBitset()];
    operations.push_back(OptimizeBranchLength{
      GetPlvIndex(P, dag_nodes_.size(), child_id),
      GetPlvIndex(R, dag_nodes_.size(), node_id),
      pcsp_idx
    });
    // Update p_hat(s)
    operations.push_back(WeightedSumAccumulate{
      GetPlvIndex(P_HAT, dag_nodes_.size(), node_id),
      pcsp_idx,
      GetPlvIndex(P, dag_nodes_.size(), child_id),
    });
  }
  // Update r_tilde(t) = r_hat(t) \circ p_hat(t).
  operations.push_back(Multiply{
    GetPlvIndex(R_TILDE, dag_nodes_.size(), node_id),
    GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
    GetPlvIndex(P_HAT, dag_nodes_.size(), node_id)
  });

  operations.push_back(Zero{GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleBranchLengthOptimization(child_id, visited_nodes, operations);
    }
    auto child_node = dag_nodes_[child_id];
    size_t pcsp_idx = pcsp_indexer_[node->GetBitset().RotateSubsplit() + child_node->GetBitset()];
    operations.push_back(OptimizeBranchLength{
      GetPlvIndex(P, dag_nodes_.size(), child_id),
      GetPlvIndex(R_TILDE, dag_nodes_.size(), node_id),
      pcsp_idx
    });
    operations.push_back(WeightedSumAccumulate{
      GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id),
      pcsp_idx,
      GetPlvIndex(P, dag_nodes_.size(), child_id),
    });
  }
  // Update r(t) = r_hat(t) \circ p_hat_tilde(t).
  operations.push_back(Multiply{
    GetPlvIndex(R, dag_nodes_.size(), node_id),
    GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
    GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id)
  });

  // Update p(t).
  operations.push_back(Multiply{
    GetPlvIndex(P, dag_nodes_.size(), node_id),
    GetPlvIndex(P_HAT, dag_nodes_.size(), node_id),
    GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id)
  });
};

GPOperationVector GPInstance::BranchLengthOptimization()
{
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (auto rootsplit : rootsplits_) {
    ScheduleBranchLengthOptimization(subsplit_to_index_[rootsplit + ~rootsplit],
                                     visited_nodes,
                                     operations);
  }
  return operations;
}

void GPInstance::ScheduleSBNParametersOptimization(size_t node_id,
                                                   std::unordered_set<size_t> &visited_nodes,
                                                   GPOperationVector &operations) {
  visited_nodes.insert(node_id);
  auto node = dag_nodes_[node_id];

  if (!node->IsRoot()) {
    // We compute R_HAT(s) = \sum_{t : s < t} q(s|t) P(s|t) r(t).
    // This is necessary to reflect changes to r(t) as well as new values
    // for q(s|t).
    // TODO: Code refactoring -- similar code is repeated twice.
    operations.push_back(Zero{GetPlvIndex(R_HAT, dag_nodes_.size(), node_id)});
    for (size_t parent_id : node->GetRootwardSorted()) {
      auto parent_node = dag_nodes_[parent_id];
      size_t pcsp_idx = pcsp_indexer_[parent_node->GetBitset() + node->GetBitset()];
      operations.push_back(WeightedSumAccumulate{
        GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
        pcsp_idx,
        GetPlvIndex(R, dag_nodes_.size(), parent_id)
      });
    }
    for (size_t parent_id : node->GetRootwardRotated()) {
      auto parent_node = dag_nodes_[parent_id];
      size_t pcsp_idx = pcsp_indexer_[parent_node->GetBitset().RotateSubsplit() + node->GetBitset()];
      operations.push_back(WeightedSumAccumulate{
        GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
        pcsp_idx,
        GetPlvIndex(R_TILDE, dag_nodes_.size(), parent_id)
      });
    }

    // Update r(s) and r_tilde(s) using rhat(s).
    operations.push_back(Multiply{
      GetPlvIndex(R, dag_nodes_.size(), node_id),
      GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
      GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id)
    });
    operations.push_back(Multiply{
      GetPlvIndex(R_TILDE, dag_nodes_.size(), node_id),
      GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
      GetPlvIndex(P_HAT, dag_nodes_.size(), node_id)
    });
  }

  if (node->IsLeaf())
    return;

  operations.push_back(Zero{GetPlvIndex(P_HAT, dag_nodes_.size(), node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleSBNParametersOptimization(child_id, visited_nodes, operations);
    }

    auto child_node = dag_nodes_[child_id];
    size_t pcsp_idx = pcsp_indexer_[node->GetBitset() + child_node->GetBitset()];
    // Update p_hat(s)
    operations.push_back(WeightedSumAccumulate{
      GetPlvIndex(P_HAT, dag_nodes_.size(), node_id),
      pcsp_idx,
      GetPlvIndex(P, dag_nodes_.size(), child_id),
    });
    operations.push_back(Likelihood{
      pcsp_idx,
      GetPlvIndex(PlvType::R, dag_nodes_.size(), node->Id()),
      GetPlvIndex(PlvType::P, dag_nodes_.size(), child_node->Id())
    });
  }
  OptimizeSBNParameters(node->GetBitset(), operations);

  // Update r_tilde(t) = r_hat(t) \circ p_hat(t).
  operations.push_back(Multiply{
    GetPlvIndex(R_TILDE, dag_nodes_.size(), node_id),
    GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
    GetPlvIndex(P_HAT, dag_nodes_.size(), node_id)
  });

  operations.push_back(Zero{GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleSBNParametersOptimization(child_id, visited_nodes, operations);
    }
    auto child_node = dag_nodes_[child_id];
    size_t pcsp_idx = pcsp_indexer_[node->GetBitset().RotateSubsplit() + child_node->GetBitset()];
    operations.push_back(WeightedSumAccumulate{
      GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id),
      pcsp_idx,
      GetPlvIndex(P, dag_nodes_.size(), child_id),
    });
    operations.push_back(Likelihood{
      pcsp_idx,
      GetPlvIndex(PlvType::R_TILDE, dag_nodes_.size(), node->Id()),
      GetPlvIndex(PlvType::P, dag_nodes_.size(), child_node->Id())
    });
  }
  OptimizeSBNParameters(node->GetBitset().RotateSubsplit(), operations);

  // Update r(t) = r_hat(t) \circ p_hat_tilde(t).
  operations.push_back(Multiply{
    GetPlvIndex(R, dag_nodes_.size(), node_id),
    GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
    GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id)
  });

  // Update p(t).
  operations.push_back(Multiply{
    GetPlvIndex(P, dag_nodes_.size(), node_id),
    GetPlvIndex(P_HAT, dag_nodes_.size(), node_id),
    GetPlvIndex(P_HAT_TILDE, dag_nodes_.size(), node_id)
  });
};

GPOperationVector GPInstance::SBNParameterOptimization()
{
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (size_t i = 0; i < rootsplits_.size(); i++) {
    auto rootsplit = rootsplits_[i];
    ScheduleSBNParametersOptimization(subsplit_to_index_[rootsplit + ~rootsplit],
                                      visited_nodes,
                                      operations);
    auto node_id = subsplit_to_index_[rootsplit + ~rootsplit];
    operations.push_back(MarginalLikelihood{
      GetPlvIndex(R_HAT, dag_nodes_.size(), node_id),
      i,
      GetPlvIndex(P, dag_nodes_.size(), node_id)
    });
  }
  // Optimize SBN parameters for the rootsplits: at this point, p-vectors are
  // already updated in the call to ScheduleSBNParametersOptimization.
  operations.push_back(UpdateSBNProbabilities{0, rootsplits_.size()});
  return operations;
}

void GPInstance::PopulatePLVs()
{
  // Performs rootward and leafward pass to populate all relevant PLVs.
  RootwardPass(rootward_order_);
  LeafwardPass(leafward_order_);
}

void GPInstance::ComputeLikelihoods()
{
  size_t taxon_count = alignment_.SequenceCount();
  // Compute likelihood values l(s|t) for each PCSP s|t.
  GPOperationVector operations;
  for (size_t i = taxon_count; i < dag_nodes_.size(); i++) {
    auto node = dag_nodes_[i];
    for (size_t child_idx : node->GetLeafwardSorted()) {
      auto child_node = dag_nodes_[child_idx];
      auto pcsp_idx = pcsp_indexer_[node->GetBitset() + child_node->GetBitset()];
      operations.push_back(Likelihood{
        pcsp_idx,
        GetPlvIndex(PlvType::R, dag_nodes_.size(), node->Id()),
        GetPlvIndex(PlvType::P, dag_nodes_.size(), child_node->Id())
      });
    }
    for (size_t child_idx : node->GetLeafwardRotated()) {
      auto child_node = dag_nodes_[child_idx];
      auto pcsp_idx = pcsp_indexer_[node->GetBitset().RotateSubsplit() + child_node->GetBitset()];
      operations.push_back(Likelihood{
        pcsp_idx,
        GetPlvIndex(PlvType::R_TILDE, dag_nodes_.size(), node->Id()),
        GetPlvIndex(PlvType::P, dag_nodes_.size(), child_node->Id())
      });
    }
  }

  // Compute marginal likelihood.
  for (size_t root_idx = 0; root_idx < rootsplits_.size(); root_idx++) {
    auto rootsplit = rootsplits_[root_idx];
    auto node_idx = subsplit_to_index_[rootsplit + ~rootsplit];
    operations.push_back(MarginalLikelihood{
      GetPlvIndex(PlvType::R_HAT, dag_nodes_.size(), node_idx),
      root_idx,
      GetPlvIndex(PlvType::P, dag_nodes_.size(), node_idx)
    });
  }

  GetEngine()->ProcessOperations(operations);
}

GPOperationVector GPInstance::MarginalLikelihoodOperations() {
  // Compute marginal likelihood.
  GPOperationVector operations;
  for (size_t i = 0; i < rootsplits_.size(); i++) {
    auto rootsplit = rootsplits_[i];
    auto root_subsplit = rootsplit + ~rootsplit;
    size_t root_idx = subsplit_to_index_[root_subsplit];
    operations.push_back(MarginalLikelihood{
      GetPlvIndex(PlvType::R_HAT, dag_nodes_.size(), root_idx),
      i,
      GetPlvIndex(PlvType::P, dag_nodes_.size(), root_idx)
    });
  }
  return operations;
}

void GPInstance::EstimateBranchLengths(double tol, size_t max_iter) {
  std::cout << "Begin branch optimization\n";
  GPOperationVector branch_optimization_operations = BranchLengthOptimization();
  GPOperationVector marginal_lik_operations = MarginalLikelihoodOperations();

  SetRootwardZero();
  SetLeafwardZero();
  GetEngine()->ResetLogMarginalLikelihood();
  PopulatePLVs();
  ComputeLikelihoods();
  double current_marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();

  for (size_t i = 0; i < max_iter; i++) {
    std::cout << "Iteration: " << (i + 1) << std::endl;
    GetEngine()->ProcessOperations(branch_optimization_operations);
    GetEngine()->ResetLogMarginalLikelihood();
    GetEngine()->ProcessOperations(marginal_lik_operations);
    double marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
    std::cout << "Current marginal log likelihood: ";
    std::cout << std::setprecision(9) << current_marginal_log_lik << std::endl;
    std::cout << "New marginal log likelihood: ";
    std::cout << std::setprecision(9) << marginal_log_lik << std::endl;
    if (marginal_log_lik < current_marginal_log_lik) {
      std::cout << "Marginal log likelihood decreased. Check branch optimization routine.\n";
      //break;
    }
    if (abs(current_marginal_log_lik - marginal_log_lik) < tol) {
      std::cout << "Converged.\n";
      break;
    }
    current_marginal_log_lik = marginal_log_lik;
  }
}

void GPInstance::EstimateSBNParameters(double tol, size_t max_iter) {
  std::cout << "Begin SBN parameter optimization\n";
  GPOperationVector sbn_param_optimization_operations = SBNParameterOptimization();
  GPOperationVector marginal_lik_operations = MarginalLikelihoodOperations();

  SetRootwardZero();
  SetLeafwardZero();
  GetEngine()->ResetLogMarginalLikelihood();
  PopulatePLVs();
  ComputeLikelihoods();

  GetEngine()->ProcessOperations(sbn_param_optimization_operations);
  GetEngine()->ProcessOperations(marginal_lik_operations);
  double marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
  std::cout << std::setprecision(9) << marginal_log_lik << std::endl;
}
