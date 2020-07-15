// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_dag.hpp"
#include "numerical_utils.hpp"

using namespace GPOperations;

GPDAG::GPDAG() : taxon_count_(0), rootsplit_and_pcsp_count_(0) {}

GPDAG::GPDAG(const RootedTreeCollection &tree_collection) {
  ProcessTrees(tree_collection);
  BuildNodes();
  BuildEdges();
  BuildPCSPIndexer();
}

size_t GPDAG::NodeCount() const { return dag_nodes_.size(); }

size_t GPDAG::RootsplitAndPCSPCount() const { return rootsplit_and_pcsp_count_; }

size_t GPDAG::GeneralizedPCSPCount() const {
  size_t fake_subsplit_parameter_count = 0;
  for (size_t taxon_idx = 0; taxon_idx < taxon_count_; taxon_idx++) {
    fake_subsplit_parameter_count += dag_nodes_[taxon_idx]->GetRootwardRotated().size();
    fake_subsplit_parameter_count += dag_nodes_[taxon_idx]->GetRootwardSorted().size();
  }
  return RootsplitAndPCSPCount() + fake_subsplit_parameter_count;
}

void GPDAG::Print() const {
  for (size_t i = 0; i < dag_nodes_.size(); i++) {
    std::cout << dag_nodes_[i]->ToString() << std::endl;
  }
}

void GPDAG::PrintGPCSPIndexer() const {
  for (const auto &[pcsp, idx] : gpcsp_indexer_) {
    std::cout << pcsp.SubsplitToString() << ", " << idx << std::endl;
  }
}

GPDAGNode *GPDAG::GetDagNode(const size_t node_id) const {
  return dag_nodes_.at(node_id).get();
}

size_t GPDAG::GetPLVIndexStatic(PLVType plv_type, size_t node_count, size_t src_idx) {
  switch (plv_type) {
    case PLVType::P:
      return src_idx;
    case PLVType::P_HAT:
      return node_count + src_idx;
    case PLVType::P_HAT_TILDE:
      return 2 * node_count + src_idx;
    case PLVType::R_HAT:
      return 3 * node_count + src_idx;
    case PLVType::R:
      return 4 * node_count + src_idx;
    case PLVType::R_TILDE:
      return 5 * node_count + src_idx;
    default:
      Failwith("Invalid PLV index requested.");
  }
}

size_t GPDAG::GetPLVIndex(PLVType plv_type, size_t src_idx) const {
  return GetPLVIndexStatic(plv_type, dag_nodes_.size(), src_idx);
}

EigenVectorXd GPDAG::BuildUniformQ() const {
  EigenVectorXd q = EigenVectorXd::Ones(GeneralizedPCSPCount());
  q.segment(0, rootsplits_.size()).array() = 1. / rootsplits_.size();
  for (const auto &[_, range] : subsplit_to_range_) {
    const auto num_child_subsplits = range.second - range.first;
    double val = 1. / num_child_subsplits;
    q.segment(range.first, num_child_subsplits).array() = val;
  }
  return q;
}

GPOperationVector GPDAG::BranchLengthOptimization() const {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (auto rootsplit : rootsplits_) {
    ScheduleBranchLengthOptimization(subsplit_to_id_.at(rootsplit + ~rootsplit),
                                     visited_nodes, operations);
  }
  return operations;
}

GPOperationVector GPDAG::ComputeLikelihoods() const {
  GPOperationVector operations;
  IterateOverRealNodes([this, &operations](const GPDAGNode *node) {
    for (size_t child_idx : node->GetLeafwardSorted()) {
      const auto child_node = GetDagNode(child_idx);
      const auto gpcsp_idx =
          gpcsp_indexer_.at(node->GetBitset() + child_node->GetBitset());
      operations.push_back(Likelihood{gpcsp_idx, GetPLVIndex(PLVType::R, node->Id()),
                                      GetPLVIndex(PLVType::P, child_node->Id())});
    }
    for (size_t child_idx : node->GetLeafwardRotated()) {
      const auto child_node = GetDagNode(child_idx);
      const auto gpcsp_idx = gpcsp_indexer_.at(node->GetBitset().RotateSubsplit() +
                                               child_node->GetBitset());
      operations.push_back(Likelihood{gpcsp_idx,
                                      GetPLVIndex(PLVType::R_TILDE, node->Id()),
                                      GetPLVIndex(PLVType::P, child_node->Id())});
    }
  });

  const auto marginal_likelihood_operations = MarginalLikelihood();
  operations.insert(operations.end(), marginal_likelihood_operations.begin(),
                    marginal_likelihood_operations.end());
  return operations;
}

GPOperationVector GPDAG::LeafwardPass() const {
  return LeafwardPass(LeafwardPassTraversal());
}

GPOperationVector GPDAG::MarginalLikelihood() const {
  GPOperationVector operations;
  for (size_t rootsplit_idx = 0; rootsplit_idx < rootsplits_.size(); rootsplit_idx++) {
    const auto rootsplit = rootsplits_[rootsplit_idx];
    const auto root_subsplit = rootsplit + ~rootsplit;
    size_t root_idx = subsplit_to_id_.at(root_subsplit);
    operations.push_back(GPOperations::IncrementMarginalLikelihood{
        GetPLVIndex(PLVType::R_HAT, root_idx), rootsplit_idx,
        GetPLVIndex(PLVType::P, root_idx)});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass() const {
  return RootwardPass(RootwardPassTraversal());
}

GPOperationVector GPDAG::OptimizeSBNParameters() const {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (size_t &id : LeafwardPassTraversal()) {
    const auto node = GetDagNode(id);
    OptimizeSBNParametersForASubsplit(node->GetBitset(), operations);
    OptimizeSBNParametersForASubsplit(node->GetBitset().RotateSubsplit(), operations);
  }
  operations.push_back(UpdateSBNProbabilities{0, rootsplits_.size()});
  return operations;
}

GPOperationVector GPDAG::SetLeafwardZero() const {
  GPOperationVector operations;
  const auto node_count = dag_nodes_.size();
  for (size_t i = 0; i < node_count; i++) {
    operations.push_back(Zero{GetPLVIndex(PLVType::R_HAT, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::R, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::R_TILDE, i)});
  }
  for (const auto &rootsplit : rootsplits_) {
    size_t root_id = subsplit_to_id_.at(rootsplit + ~rootsplit);
    operations.push_back(
        SetToStationaryDistribution{GetPLVIndex(PLVType::R_HAT, root_id)});
  }
  return operations;
}

GPOperationVector GPDAG::SetRhatToStationary() const {
  GPOperationVector operations;
  for (const auto &rootsplit : rootsplits_) {
    size_t root_id = subsplit_to_id_.at(rootsplit + ~rootsplit);
    operations.push_back(
        SetToStationaryDistribution{GetPLVIndex(PLVType::R_HAT, root_id)});
  }
  return operations;
}

GPOperationVector GPDAG::SetRootwardZero() const {
  GPOperationVector operations;
  const auto node_count = dag_nodes_.size();
  for (size_t i = taxon_count_; i < node_count; i++) {
    operations.push_back(Zero{GetPLVIndex(PLVType::P, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT_TILDE, i)});
  }
  return operations;
}

void GPDAG::IterateOverRealNodes(std::function<void(const GPDAGNode *)> f) const {
  Assert(taxon_count_ < dag_nodes_.size(), "No real DAG nodes!");
  for (auto it = dag_nodes_.cbegin() + taxon_count_; it < dag_nodes_.cend(); it++) {
    f((*it).get());
  }
}

std::vector<Bitset> GPDAG::GetChildrenSubsplits(const Bitset &subsplit,
                                                bool include_fake_subsplits) {
  std::vector<Bitset> children_subsplits;

  if (parent_to_range_.count(subsplit)) {
    const auto [start, stop] = parent_to_range_.at(subsplit);
    for (auto idx = start; idx < stop; idx++) {
      const auto child_subsplit = index_to_child_.at(idx);
      children_subsplits.push_back(child_subsplit);
    }
  } else if (include_fake_subsplits) {
    // In the case where second chunk of the subsplit is just a single taxon,
    // the subsplit will not map to any value (parent_to_range_[subsplit] doesn't
    // exist). But we still need to create and connect to fake subsplits in the DAG. A
    // subsplit has a fake subsplit as a child if the first chunk is non-zero and the
    // second chunk has exactly one bit set to 1.
    if (subsplit.SplitChunk(0).Any() &&
        subsplit.SplitChunk(1).SingletonOption().has_value()) {
      // The fake subsplit corresponds to the second chunk of subsplit.
      // Prepend it by 0's.
      Bitset zero(subsplit.size() / 2);
      Bitset fake_subsplit = zero + subsplit.SplitChunk(1);
      children_subsplits.push_back(fake_subsplit);
    }
  }

  return children_subsplits;
}

void GPDAG::ProcessTrees(const RootedTreeCollection &tree_collection) {
  taxon_count_ = tree_collection.TaxonCount();
  const auto topology_counter = tree_collection.TopologyCounter();

  std::tie(rootsplits_, std::ignore, index_to_child_, parent_to_range_,
           rootsplit_and_pcsp_count_) =
      SBNMaps::BuildIndexerBundle(RootedSBNMaps::RootsplitCounterOf(topology_counter),
                                  RootedSBNMaps::PCSSCounterOf(topology_counter));
}

void GPDAG::CreateAndInsertNode(const Bitset &subsplit) {
  size_t id = dag_nodes_.size();
  SafeInsert(subsplit_to_id_, subsplit, id);
  dag_nodes_.push_back(std::make_unique<GPDAGNode>(id, subsplit));
}

void GPDAG::ConnectNodes(size_t idx, bool rotated) {
  const auto node = dag_nodes_[idx].get();
  // Retrieve children subsplits, set edge relation.
  const Bitset subsplit = node->GetBitset(rotated);
  const auto children = GetChildrenSubsplits(subsplit, true);
  for (const auto &child_subsplit : children) {
    const auto child_node = GetDagNode(subsplit_to_id_.at(child_subsplit));
    if (rotated) {
      node->AddLeafwardRotated(child_node->Id());
      child_node->AddRootwardRotated(node->Id());
    } else {
      node->AddLeafwardSorted(child_node->Id());
      child_node->AddRootwardSorted(node->Id());
    }
  }
}

void GPDAG::BuildNodesDepthFirst(const Bitset &subsplit,
                                 std::unordered_set<Bitset> &visited_subsplits) {
  visited_subsplits.insert(subsplit);
  auto children_subsplits = GetChildrenSubsplits(subsplit, false);
  for (const auto &child_subsplit : children_subsplits) {
    if (!visited_subsplits.count(child_subsplit)) {
      BuildNodesDepthFirst(child_subsplit, visited_subsplits);
    }
  }
  children_subsplits = GetChildrenSubsplits(subsplit.RotateSubsplit(), false);
  for (const auto &child_subsplit : children_subsplits) {
    if (!visited_subsplits.count(child_subsplit)) {
      BuildNodesDepthFirst(child_subsplit, visited_subsplits);
    }
  }
  CreateAndInsertNode(subsplit);
}

void GPDAG::BuildNodes() {
  std::unordered_set<Bitset> visited_subsplits;

  // We will create fake subsplits and insert to dag_nodes_.
  // These nodes will take IDs in [0, taxon_count_).
  Bitset zero(taxon_count_);
  for (size_t taxon_idx = 0; taxon_idx < taxon_count_; taxon_idx++) {
    Bitset fake(taxon_count_);
    fake.set(taxon_idx);
    Bitset fake_subsplit = zero + fake;
    CreateAndInsertNode(fake_subsplit);
  }
  // We are going to add the remaining nodes.
  // The root splits will take on the higher IDs compared to the non-rootsplits.
  for (const auto &rootsplit : rootsplits_) {
    const auto subsplit = rootsplit + ~rootsplit;
    BuildNodesDepthFirst(subsplit, visited_subsplits);
  }
}

void GPDAG::BuildEdges() {
  for (size_t i = taxon_count_; i < dag_nodes_.size(); i++) {
    ConnectNodes(i, false);
    ConnectNodes(i, true);
  }
}

void GPDAG::BuildPCSPIndexer() {
  size_t idx = 0;
  for (const auto &rootsplit : rootsplits_) {
    SafeInsert(gpcsp_indexer_, rootsplit + ~rootsplit, idx);
    idx++;
  }

  IterateOverRealNodes([this, &idx](const GPDAGNode *node) {
    auto child_count = node->GetLeafwardSorted().size();
    if (child_count > 0) {
      SafeInsert(subsplit_to_range_, node->GetBitset(), {idx, idx + child_count});
      for (size_t j = 0; j < child_count; j++) {
        const auto child = GetDagNode(node->GetLeafwardSorted().at(j));
        SafeInsert(gpcsp_indexer_, node->GetBitset() + child->GetBitset(), idx);
        idx++;
      }
    }
    child_count = node->GetLeafwardRotated().size();
    if (child_count > 0) {
      SafeInsert(subsplit_to_range_, node->GetBitset().RotateSubsplit(),
                 {idx, idx + child_count});
      for (size_t j = 0; j < child_count; j++) {
        const auto child = GetDagNode(node->GetLeafwardRotated().at(j));
        SafeInsert(gpcsp_indexer_,
                   node->GetBitset().RotateSubsplit() + child->GetBitset(), idx);
        idx++;
      }
    }
  });
}

void RootwardDepthFirst(size_t id,
                        const std::vector<std::unique_ptr<GPDAGNode>> &dag_nodes,
                        std::vector<size_t> &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  SafeInsert(visited_nodes, id);
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
                        const std::vector<std::unique_ptr<GPDAGNode>> &dag_nodes,
                        std::vector<size_t> &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  SafeInsert(visited_nodes, id);
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

GPOperationVector GPDAG::LeafwardPass(std::vector<size_t> visit_order) const {
  GPOperationVector operations;
  for (const size_t node_id : visit_order) {
    const auto node = GetDagNode(node_id);

    // Build rhat(s) via rhat(s) += \sum_t q(s|t) P'(s|t) r(t)
    AddRhatOperations(node, operations);
    // Multiply to get r(s) = rhat(s) \circ phat(s_tilde).
    operations.push_back(Multiply{GetPLVIndex(PLVType::R, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
    // Multiply to get r(s_tilde) = rhat(s) \circ phat(s).
    operations.push_back(Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT, node_id)});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass(std::vector<size_t> visit_order) const {
  GPOperationVector operations;
  for (const size_t node_id : visit_order) {
    const auto node = GetDagNode(node_id);
    if (node->IsLeaf()) {
      continue;
    }
    // Build phat(s).
    AddPhatOperations(node, false, operations);
    // Build phat(s_tilde).
    AddPhatOperations(node, true, operations);
    // Multiply to get p(s) = phat(s) \circ phat(s_tilde).
    operations.push_back(Multiply{node_id, GetPLVIndex(PLVType::P_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
  }
  return operations;
}

std::vector<size_t> GPDAG::LeafwardPassTraversal() const {
  std::vector<size_t> visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (size_t leaf_idx = 0; leaf_idx < taxon_count_; leaf_idx++) {
    RootwardDepthFirst(leaf_idx, dag_nodes_, visit_order, visited_nodes);
  }
  return visit_order;
}

std::vector<size_t> GPDAG::RootwardPassTraversal() const {
  std::vector<size_t> visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (const auto &rootsplit : rootsplits_) {
    size_t root_idx = subsplit_to_id_.at(rootsplit + ~rootsplit);
    LeafwardDepthFirst(root_idx, dag_nodes_, visit_order, visited_nodes);
  }
  return visit_order;
}

void GPDAG::AddPhatOperations(const GPDAGNode *node, bool rotated,
                              GPOperationVector &operations) const {
  std::vector<size_t> child_idxs =
      rotated ? node->GetLeafwardRotated() : node->GetLeafwardSorted();
  PLVType plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
  const auto parent_subsplit = node->GetBitset(rotated);
  for (size_t child_idx : child_idxs) {
    const auto child_node = GetDagNode(child_idx);
    const auto child_subsplit = child_node->GetBitset();
    const auto pcsp = parent_subsplit + child_subsplit;
    if (!gpcsp_indexer_.count(pcsp)) {
      Failwith("Non-existent PCSP index.");
    }
    const auto gpcsp_idx = gpcsp_indexer_.at(pcsp);

    operations.push_back(
        IncrementWithWeightedEvolvedPLV{GetPLVIndex(plv_type, node->Id()), gpcsp_idx,
                                        GetPLVIndex(PLVType::P, child_idx)});
  }
}

void GPDAG::AddRhatOperations(const GPDAGNode *node,
                              GPOperationVector &operations) const {
  const auto subsplit = node->GetBitset();
  for (const size_t parent_idx : node->GetRootwardSorted()) {
    const auto parent_node = GetDagNode(parent_idx);
    const auto parent_subsplit = parent_node->GetBitset();
    const auto gpcsp_idx = gpcsp_indexer_.at(parent_subsplit + subsplit);

    operations.push_back(IncrementWithWeightedEvolvedPLV{
        GetPLVIndex(PLVType::R_HAT, node->Id()), gpcsp_idx,
        GetPLVIndex(PLVType::R, parent_node->Id())});
  }
  for (const size_t parent_idx : node->GetRootwardRotated()) {
    const auto parent_node = GetDagNode(parent_idx);
    const auto parent_subsplit = parent_node->GetBitset().RotateSubsplit();
    const auto gpcsp_idx = gpcsp_indexer_.at(parent_subsplit + subsplit);

    operations.push_back(IncrementWithWeightedEvolvedPLV{
        GetPLVIndex(PLVType::R_HAT, node->Id()), gpcsp_idx,
        GetPLVIndex(PLVType::R_TILDE, parent_node->Id())});
  }
}

void GPDAG::OptimizeSBNParametersForASubsplit(const Bitset &subsplit,
                                              GPOperationVector &operations) const {
  if (subsplit_to_range_.count(subsplit)) {
    const auto param_range = subsplit_to_range_.at(subsplit);
    if (param_range.second - param_range.first > 1) {
      operations.push_back(
          UpdateSBNProbabilities{param_range.first, param_range.second});
    }
  }
}

void GPDAG::ScheduleBranchLengthOptimization(size_t node_id,
                                             std::unordered_set<size_t> &visited_nodes,
                                             GPOperationVector &operations) const {
  visited_nodes.insert(node_id);
  const auto node = GetDagNode(node_id);

  if (!node->IsRoot()) {
    // Compute R_HAT(s) = \sum_{t : s < t} q(s|t) P(s|t) r(t).
    operations.push_back(Zero{GetPLVIndex(PLVType::R_HAT, node_id)});
    UpdateRHat(node_id, false, operations);
    UpdateRHat(node_id, true, operations);

    // Update r(s) and r_tilde(s)
    operations.push_back(Multiply{GetPLVIndex(PLVType::R, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
    operations.push_back(Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT, node_id)});
  }

  if (node->IsLeaf()) return;

  operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT, node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleBranchLengthOptimization(child_id, visited_nodes, operations);
    }
    OptimizeBranchLengthUpdatePHat(node_id, child_id, false, operations);
  }
  // Update r_tilde(t) = r_hat(t) \circ p_hat(t).
  operations.push_back(Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                                GetPLVIndex(PLVType::R_HAT, node_id),
                                GetPLVIndex(PLVType::P_HAT, node_id)});

  operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleBranchLengthOptimization(child_id, visited_nodes, operations);
    }
    OptimizeBranchLengthUpdatePHat(node_id, child_id, true, operations);
  }
  // Update r(t) = r_hat(t) \circ p_hat_tilde(t).
  operations.push_back(Multiply{GetPLVIndex(PLVType::R, node_id),
                                GetPLVIndex(PLVType::R_HAT, node_id),
                                GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});

  // Update p(t).
  operations.push_back(Multiply{GetPLVIndex(PLVType::P, node_id),
                                GetPLVIndex(PLVType::P_HAT, node_id),
                                GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
};

void GPDAG::UpdateRHat(size_t node_id, bool rotated,
                       GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  PLVType src_plv_type = rotated ? PLVType::R_TILDE : PLVType::R;
  const auto parent_nodes =
      rotated ? node->GetRootwardRotated() : node->GetRootwardSorted();
  for (size_t parent_id : parent_nodes) {
    const auto parent_node = GetDagNode(parent_id);
    auto pcsp =
        rotated ? parent_node->GetBitset().RotateSubsplit() : parent_node->GetBitset();
    pcsp = pcsp + node->GetBitset();
    size_t gpcsp_idx = gpcsp_indexer_.at(pcsp);
    operations.push_back(
        IncrementWithWeightedEvolvedPLV{GetPLVIndex(PLVType::R_HAT, node_id), gpcsp_idx,
                                        GetPLVIndex(src_plv_type, parent_id)});
  }
}

void GPDAG::UpdatePHatComputeLikelihood(size_t node_id, size_t child_node_id,
                                        bool rotated,
                                        GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  const auto child_node = GetDagNode(child_node_id);
  auto pcsp = rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  pcsp = pcsp + child_node->GetBitset();
  size_t gpcsp_idx = gpcsp_indexer_.at(pcsp);
  // Update p_hat(s)
  operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  operations.push_back(Likelihood{
      gpcsp_idx, GetPLVIndex(rotated ? PLVType::R_TILDE : PLVType::R, node->Id()),
      GetPLVIndex(PLVType::P, child_node->Id())});
}

void GPDAG::OptimizeBranchLengthUpdatePHat(size_t node_id, size_t child_node_id,
                                           bool rotated,
                                           GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  const auto child_node = GetDagNode(child_node_id);
  auto pcsp = rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  pcsp = pcsp + child_node->GetBitset();
  size_t gpcsp_idx = gpcsp_indexer_.at(pcsp);
  operations.push_back(OptimizeBranchLength{
      GetPLVIndex(PLVType::P, child_node_id),
      GetPLVIndex(rotated ? PLVType::R_TILDE : PLVType::R, node_id), gpcsp_idx});
  // Update p_hat(s)
  operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
}
