// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_dag.hpp"

#include <iostream>

#include "numerical_utils.hpp"

using namespace GPOperations;  // NOLINT

GPDAG::GPDAG()
    : taxon_count_(0), gpcsp_count_without_fake_subsplits_(0), topology_count_(0.) {}

GPDAG::GPDAG(const RootedTreeCollection &tree_collection) : GPDAG() {
  ProcessTrees(tree_collection);
  BuildNodes();
  BuildEdges();
  ExpandRootsplitsInIndexer();
  AddFakeSubsplitsToGPCSPIndexerAndParentToRange();
  CountTopologies();
}

void GPDAG::CountTopologies() {
  topology_count_below_ = EigenVectorXd::Ones(NodeCount());
  for (const auto &node_id : RootwardPassTraversal()) {
    const auto &node = GetDagNode(node_id);
    if (!node->IsLeaf()) {
      for (const bool rotated : {false, true}) {
        double per_rotated_count = 0.;
        // Sum options across the possible children.
        for (const auto &child_id : node->GetLeafward(rotated)) {
          per_rotated_count += topology_count_below_[child_id];
        }
        // Take the product across the number of options for the left and right branches
        // of the tree.
        topology_count_below_[node_id] *= per_rotated_count;
      }
    }
  }
  topology_count_ = 0;
  IterateOverRootsplitIds(
      [this](size_t root_id) { topology_count_ += topology_count_below_[root_id]; });
}

double GPDAG::TopologyCount() const { return topology_count_; }

size_t GPDAG::NodeCount() const { return dag_nodes_.size(); }

size_t GPDAG::RootsplitCount() const { return rootsplits_.size(); }
size_t GPDAG::GPCSPCount() const { return gpcsp_count_without_fake_subsplits_; }
size_t GPDAG::GPCSPCountWithFakeSubsplits() const { return gpcsp_indexer_.size(); }

void GPDAG::Print() const {
  for (const auto &dag_node : dag_nodes_) {
    std::cout << dag_node->ToString() << std::endl;
  }
}

void GPDAG::PrintGPCSPIndexer() const {
  for (const auto &[pcsp, idx] : gpcsp_indexer_) {
    // Since pcsp may be a rootsubsplit, we need to check and get subsplit.
    std::string gpcsp_string = pcsp.size() == (taxon_count_ * 3)
                                   ? pcsp.PCSPToString()
                                   : pcsp.SubsplitToString();
    std::cout << gpcsp_string << ", " << idx << std::endl;
  }
}

const BitsetSizeMap &GPDAG::GetGPCSPIndexer() const { return gpcsp_indexer_; }
const BitsetSizePairMap &GPDAG::GetParentToRange() const { return parent_to_range_; }

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

size_t GPDAG::GetGPCSPIndex(const Bitset &parent_subsplit,
                            const Bitset &child_subsplit) const {
  return gpcsp_indexer_.at(Bitset::PCSPOfPair(parent_subsplit, child_subsplit));
}

EigenVectorXd GPDAG::BuildUniformQ() const {
  EigenVectorXd q = EigenVectorXd::Ones(GPCSPCountWithFakeSubsplits());
  q.segment(0, rootsplits_.size()).array() = 1. / rootsplits_.size();
  for (const auto &[_, range] : parent_to_range_) {
    const auto num_child_subsplits = range.second - range.first;
    double val = 1. / num_child_subsplits;
    q.segment(range.first, num_child_subsplits).array() = val;
  }
  return q;
}

Node::NodePtrVec GPDAG::GenerateAllGPNodeIndexedTopologies() const {
  std::vector<Node::NodePtrVec> topology_below(NodeCount());

  auto GetSubtopologies = [&topology_below](const GPDAGNode *node) {
    Node::NodePtrVec rotated_subtopologies, sorted_subtopologies;
    for (const bool rotated : {false, true}) {
      for (const auto &child_id : node->GetLeafward(rotated)) {
        for (const auto &subtopology : topology_below.at(child_id)) {
          rotated ? rotated_subtopologies.push_back(subtopology)
                  : sorted_subtopologies.push_back(subtopology);
        }
      }
    }
    return std::make_pair(rotated_subtopologies, sorted_subtopologies);
  };

  auto MergeTopologies = [](size_t node_id, Node::NodePtrVec &rotated_subtopologies,
                            Node::NodePtrVec &sorted_subtopologies) {
    Node::NodePtrVec topologies;
    for (const auto &rotated_subtopology : rotated_subtopologies) {
      for (const auto &sorted_subtopology : sorted_subtopologies) {
        Node::NodePtr new_topology =
            Node::Join(sorted_subtopology, rotated_subtopology, node_id);
        topologies.push_back(new_topology);
      }
    }
    return topologies;
  };

  for (const auto &node_id : RootwardPassTraversal()) {
    const auto &node = GetDagNode(node_id);
    if (node->IsLeaf()) {
      topology_below.at(node_id).push_back(Node::Leaf(node_id));
    } else {
      auto [rotated_topologies, sorted_topologies] = GetSubtopologies(node);
      topology_below[node_id] =
          MergeTopologies(node_id, rotated_topologies, sorted_topologies);
    }
  }

  Node::NodePtrVec topologies;
  IterateOverRootsplitIds([&topologies, &topology_below](size_t root_id) {
    topologies.insert(topologies.end(), topology_below.at(root_id).begin(),
                      topology_below.at(root_id).end());
  });

  Assert(topologies.size() == TopologyCount(),
         "The realized number of topologies does not match the expected count.");

  return topologies;
}

EigenVectorXd GPDAG::BuildUniformPrior() const {
  EigenVectorXd q = EigenVectorXd::Ones(GPCSPCountWithFakeSubsplits());

  for (const auto &node_id : RootwardPassTraversal()) {
    const auto &node = GetDagNode(node_id);
    if (!node->IsLeaf()) {
      for (const bool rotated : {false, true}) {
        double per_rotated_count = 0.;
        for (const auto &child_id : node->GetLeafward(rotated)) {
          per_rotated_count += topology_count_below_[child_id];
        }
        for (const auto &child_id : node->GetLeafward(rotated)) {
          size_t gpcsp_idx = GetGPCSPIndex(node->GetBitset(rotated),
                                           GetDagNode(child_id)->GetBitset());
          q(gpcsp_idx) = topology_count_below_(child_id) / per_rotated_count;
        }
      }
    }
  }

  IterateOverRootsplitIds([this, &q](size_t root_id) {
    auto node = GetDagNode(root_id);
    auto gpcsp_idx = gpcsp_indexer_.at(node->GetBitset());
    q[gpcsp_idx] = topology_count_below_[root_id] / topology_count_;
  });

  return q;
}

void GPDAG::UpdateRPLVs(size_t node_id, GPOperationVector &operations) const {
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

void GPDAG::OptimizeBranchLengthsUpdatePHatAndPropagateRPLV(
    const GPDAGNode *node, bool rotated, GPOperationVector &operations) const {
  PLVType p_hat_plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
  PLVType r_plv_type = rotated ? PLVType::R : PLVType::R_TILDE;

  size_t node_id = node->Id();

  operations.push_back(Zero{GetPLVIndex(p_hat_plv_type, node_id)});
  IterateOverLeafwardEdges(
      node, rotated,
      [this, &operations, &rotated, &node_id](const GPDAGNode *child_node) {
        OptimizeBranchLengthUpdatePHat(node_id, child_node->Id(), rotated, operations);
      });
  // Update r_tilde(t) = r_hat(t) \circ p_hat(t) and r(t) = r_hat(t) \circ
  // p_hat_tilde(t)
  operations.push_back(Multiply{GetPLVIndex(r_plv_type, node_id),
                                GetPLVIndex(PLVType::R_HAT, node_id),
                                GetPLVIndex(p_hat_plv_type, node_id)});
}

void GPDAG::ScheduleBranchLengthOptimization(bool is_reverse_postorder,
                                             GPOperationVector &operations) const {
  auto traversal = RootwardPassTraversal();
  if (is_reverse_postorder) {
    std::reverse(traversal.begin(), traversal.end());
  }

  for (auto node_id : traversal) {
    const auto node = GetDagNode(node_id);

    if (!node->IsRoot()) {
      UpdateRPLVs(node_id, operations);
    }
    if (node->IsLeaf()) {
      continue;
    }

    OptimizeBranchLengthsUpdatePHatAndPropagateRPLV(node, true, operations);
    OptimizeBranchLengthsUpdatePHatAndPropagateRPLV(node, false, operations);

    // Update p(t).
    operations.push_back(Multiply{GetPLVIndex(PLVType::P, node_id),
                                  GetPLVIndex(PLVType::P_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
  }
};

GPOperationVector GPDAG::BranchLengthOptimization() const {
  GPOperationVector operations;

  ScheduleBranchLengthOptimization(true, operations);
  ScheduleBranchLengthOptimization(false, operations);

  return operations;
}

GPOperationVector GPDAG::ComputeLikelihoods() const {
  GPOperationVector operations;
  IterateOverRealNodes([this, &operations](const GPDAGNode *node) {
    IterateOverLeafwardEdges(
        node,
        [this, node, &operations](const bool rotated, const GPDAGNode *child_node) {
          const auto gpcsp_idx =
              GetGPCSPIndex(node->GetBitset(rotated), child_node->GetBitset());
          operations.push_back(Likelihood{gpcsp_idx,
                                          GetPLVIndex(RPLVType(rotated), node->Id()),
                                          GetPLVIndex(PLVType::P, child_node->Id())});
        });
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
  return operations;
}

GPOperationVector GPDAG::SetRhatToStationary() const {
  GPOperationVector operations;
  IterateOverRootsplitIds([this, &operations](size_t rootsplit_id) {
    size_t root_gpcsp_idx = gpcsp_indexer_.at(GetDagNode(rootsplit_id)->GetBitset());
    operations.push_back(SetToStationaryDistribution{
        GetPLVIndex(PLVType::R_HAT, rootsplit_id), root_gpcsp_idx});
  });
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

void GPDAG::IterateOverRealNodes(const NodeLambda &f) const {
  Assert(taxon_count_ < dag_nodes_.size(), "No real DAG nodes!");
  for (auto it = dag_nodes_.cbegin() + taxon_count_; it < dag_nodes_.cend(); it++) {
    f((*it).get());
  }
}

void GPDAG::IterateOverLeafwardEdges(const GPDAGNode *node, bool rotated,
                                     const NodeLambda &f) const {
  for (const size_t child_idx : node->GetLeafward(rotated)) {
    f(GetDagNode(child_idx));
  }
}

void GPDAG::IterateOverLeafwardEdges(const GPDAGNode *node,
                                     const EdgeDestinationLambda &f) const {
  for (bool rotated : {false, true}) {
    for (const size_t child_idx : node->GetLeafward(rotated)) {
      f(rotated, GetDagNode(child_idx));
    }
  }
}

void GPDAG::IterateOverRootwardEdges(const GPDAGNode *node,
                                     const EdgeDestinationLambda &f) const {
  for (bool rotated : {false, true}) {
    for (const size_t parent_idx : node->GetRootward(rotated)) {
      f(rotated, GetDagNode(parent_idx));
    }
  }
}

void GPDAG::IterateOverRootsplitIds(const std::function<void(size_t)> &f) const {
  for (const auto &rootsplit : rootsplits_) {
    f(subsplit_to_id_.at(rootsplit + ~rootsplit));
  }
}

std::vector<Bitset> GPDAG::GetChildrenSubsplits(const Bitset &subsplit,
                                                bool include_fake_subsplits) {
  std::vector<Bitset> children_subsplits;
  if (parent_to_range_.count(subsplit) > 0) {
    const auto [start, stop] = parent_to_range_.at(subsplit);
    for (auto idx = start; idx < stop; idx++) {
      children_subsplits.push_back(index_to_child_.at(idx));
    }
  } else if (include_fake_subsplits) {
    // This method is designed to be called before calling
    // AddFakeSubsplitsToGPCSPIndexerAndParentToRange. In that case, if the second chunk
    // of the subsplit is just a single taxon, the subsplit will not map to any value in
    // parent_to_range_.
    //
    // But we still need to create and connect to fake subsplits in the DAG. So, here we
    // make a fake child subsplit.
    children_subsplits.push_back(Bitset::FakeChildSubsplit(subsplit));
  }

  return children_subsplits;
}

void GPDAG::ProcessTrees(const RootedTreeCollection &tree_collection) {
  taxon_count_ = tree_collection.TaxonCount();
  const auto topology_counter = tree_collection.TopologyCounter();

  std::tie(rootsplits_, gpcsp_indexer_, index_to_child_, parent_to_range_,
           gpcsp_count_without_fake_subsplits_) =
      SBNMaps::BuildIndexerBundle(RootedSBNMaps::RootsplitSupportOf(topology_counter),
                                  RootedSBNMaps::PCSPSupportOf(topology_counter));
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
  for (bool rotated : {false, true}) {
    for (const auto &child_subsplit :
         GetChildrenSubsplits(PerhapsRotateSubsplit(subsplit, rotated), false)) {
      if (visited_subsplits.count(child_subsplit) == 0) {
        BuildNodesDepthFirst(child_subsplit, visited_subsplits);
      }
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
    CreateAndInsertNode(
        Bitset::FakeSubsplit(Bitset::Singleton(taxon_count_, taxon_idx)));
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

// This should go away when addressing #273.
void GPDAG::ExpandRootsplitsInIndexer() {
  BitsetSizeMap old_indexer;
  std::swap(gpcsp_indexer_, old_indexer);
  for (const auto &[gpcsp, index] : old_indexer) {
    if (gpcsp.size() == taxon_count_) {
      SafeInsert(gpcsp_indexer_, gpcsp + ~gpcsp, index);
    } else {
      SafeInsert(gpcsp_indexer_, gpcsp, index);
    }
  }
}

void GPDAG::AddFakeSubsplitsToGPCSPIndexerAndParentToRange() {
  for (size_t i = 0; i < taxon_count_; i++) {
    const auto current_bitset = dag_nodes_[i]->GetBitset();
    IterateOverRootwardEdges(
        GetDagNode(i),
        [this, current_bitset](const bool rotated, const GPDAGNode *node) {
          SafeInsert(parent_to_range_, node->GetBitset(rotated),
                     {gpcsp_indexer_.size(), gpcsp_indexer_.size() + 1});
          SafeInsert(gpcsp_indexer_,
                     Bitset::PCSPOfPair(node->GetBitset(rotated), current_bitset),
                     gpcsp_indexer_.size());
        });
  }
}

void RootwardDepthFirst(size_t id,
                        const std::vector<std::unique_ptr<GPDAGNode>> &dag_nodes,
                        SizeVector &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  SafeInsert(visited_nodes, id);
  for (size_t child_id : dag_nodes.at(id)->GetRootwardSorted()) {
    if (visited_nodes.count(child_id) == 0) {
      RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : dag_nodes.at(id)->GetRootwardRotated()) {
    if (visited_nodes.count(child_id) == 0) {
      RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

void LeafwardDepthFirst(size_t id,
                        const std::vector<std::unique_ptr<GPDAGNode>> &dag_nodes,
                        SizeVector &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  SafeInsert(visited_nodes, id);
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardSorted()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardRotated()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

GPOperationVector GPDAG::LeafwardPass(const SizeVector &visit_order) const {
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

GPOperationVector GPDAG::RootwardPass(const SizeVector &visit_order) const {
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

SizeVector GPDAG::LeafwardPassTraversal() const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (size_t leaf_id = 0; leaf_id < taxon_count_; leaf_id++) {
    RootwardDepthFirst(leaf_id, dag_nodes_, visit_order, visited_nodes);
  }
  return visit_order;
}

SizeVector GPDAG::RootwardPassTraversal() const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  IterateOverRootsplitIds([this, &visit_order, &visited_nodes](size_t root_id) {
    LeafwardDepthFirst(root_id, dag_nodes_, visit_order, visited_nodes);
  });
  return visit_order;
}

// Take in some new operations, determine an appropriate PrepForMarginalization for
// them, then append the PrepForMarginalization and the new operations to `operations`
// (in that order).
void AppendOperationsAfterPrepForMarginalization(
    GPOperationVector &operations, const GPOperationVector &new_operations) {
  if (!new_operations.empty()) {
    operations.push_back(PrepForMarginalizationOfOperations(new_operations));
    operations.insert(operations.end(), new_operations.begin(), new_operations.end());
  }
}

void GPDAG::AddPhatOperations(const GPDAGNode *node, bool rotated,
                              GPOperationVector &operations) const {
  PLVType plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
  const auto parent_subsplit = node->GetBitset(rotated);
  const size_t dest_idx = GetPLVIndex(plv_type, node->Id());
  GPOperationVector new_operations;
  for (const size_t &child_idx : node->GetLeafward(rotated)) {
    const auto gpcsp_idx =
        GetGPCSPIndex(parent_subsplit, GetDagNode(child_idx)->GetBitset());
    new_operations.push_back(IncrementWithWeightedEvolvedPLV{
        dest_idx, gpcsp_idx, GetPLVIndex(PLVType::P, child_idx)});
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::AddRhatOperations(const GPDAGNode *node,
                              GPOperationVector &operations) const {
  const auto subsplit = node->GetBitset();
  GPOperationVector new_operations;
  IterateOverRootwardEdges(node, [this, node, &new_operations, subsplit](
                                     const bool rotated, const GPDAGNode *parent_node) {
    const auto parent_subsplit = parent_node->GetBitset(rotated);
    const auto gpcsp_idx = GetGPCSPIndex(parent_subsplit, subsplit);

    new_operations.push_back(IncrementWithWeightedEvolvedPLV{
        GetPLVIndex(PLVType::R_HAT, node->Id()), gpcsp_idx,
        GetPLVIndex(RPLVType(rotated), parent_node->Id())});
  });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeSBNParametersForASubsplit(const Bitset &subsplit,
                                              GPOperationVector &operations) const {
  if (parent_to_range_.count(subsplit) > 0) {
    const auto param_range = parent_to_range_.at(subsplit);
    if (param_range.second - param_range.first > 1) {
      operations.push_back(
          UpdateSBNProbabilities{param_range.first, param_range.second});
    }
  }
}

void GPDAG::UpdateRHat(size_t node_id, bool rotated,
                       GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  PLVType src_plv_type = rotated ? PLVType::R_TILDE : PLVType::R;
  const auto parent_nodes =
      rotated ? node->GetRootwardRotated() : node->GetRootwardSorted();
  GPOperationVector new_operations;
  for (size_t parent_id : parent_nodes) {
    const auto parent_node = GetDagNode(parent_id);
    auto parent_subsplit =
        rotated ? parent_node->GetBitset().RotateSubsplit() : parent_node->GetBitset();
    size_t gpcsp_idx = GetGPCSPIndex(parent_subsplit, node->GetBitset());
    new_operations.push_back(
        IncrementWithWeightedEvolvedPLV{GetPLVIndex(PLVType::R_HAT, node_id), gpcsp_idx,
                                        GetPLVIndex(src_plv_type, parent_id)});
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::UpdatePHatComputeLikelihood(size_t node_id, size_t child_node_id,
                                        bool rotated,
                                        GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  const auto child_node = GetDagNode(child_node_id);
  auto parent_subsplit =
      rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  size_t gpcsp_idx = GetGPCSPIndex(parent_subsplit, child_node->GetBitset());
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  new_operations.push_back(Likelihood{gpcsp_idx,
                                      GetPLVIndex(RPLVType(rotated), node->Id()),
                                      GetPLVIndex(PLVType::P, child_node->Id())});
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeBranchLengthUpdatePHat(size_t node_id, size_t child_node_id,
                                           bool rotated,
                                           GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  const auto child_node = GetDagNode(child_node_id);
  auto parent_subsplit =
      rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  size_t gpcsp_idx = GetGPCSPIndex(parent_subsplit, child_node->GetBitset());
  operations.push_back(OptimizeBranchLength{GetPLVIndex(PLVType::P, child_node_id),
                                            GetPLVIndex(RPLVType(rotated), node_id),
                                            gpcsp_idx});
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

Bitset GPDAG::PerhapsRotateSubsplit(const Bitset &subsplit, bool rotated) {
  return rotated ? subsplit.RotateSubsplit() : subsplit;
}
